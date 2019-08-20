from typing import List, Optional, Tuple
import logging
import os.path
import pandas as pd
from pandas import DataFrame
from blacksheep.classes import OutlierTable, qValues
from blacksheep._outlierTable import _convert_to_outliers
from blacksheep._outlierTable import _convert_to_counts
from blacksheep.comparisons import _compare_groups
from blacksheep.comparisons import get_sample_lists
from blacksheep._constants import *


def make_outliers_table(
    df: DataFrame,
    iqrs: float = 1.5,
    up_or_down: str = "up",
    aggregate: bool = True,
    save_outlier_table: bool = False,
    save_frac_table: bool = False,
    output_prefix: str = "outliers",
    ind_sep: str = "-",
) -> OutlierTable:
    """Converts a DataFrame of values into an OutlierTable object, which includes a DataFrame
    of outlier and non-outlier count values.

    Args:
        df: Input DataFrame with samples as columns and sites/genes as columns.
        iqrs: The number of inter-quartile ranges (IQRs) above or below the median to consider a \
        value as an outlier.
        up_or_down: Whether to call up or down outliers. Up is above the median; down \
        is below the median. Options "up" or "down".
        aggregate: Whether to sum outliers across a grouping (e.g. gene-level) than individual \
        sites. For instance if columns indicate phosphosites on proteins, with the format \
        "RAG2-S365", output will show counts of outliers per protein (e.g. RAG2) rather than on \
        individual sites (e.g. RAG2-S365).
        save_outlier_table: Whether to write a file with the outlier count table.
        save_frac_table: Whether to write a file with the outlier fraction table.
        output_prefix: If files are written, a prefix for the files.
        ind_sep: The separator used in sites, for instance, to separate a gene and site. \
        If just using genes (i.e. no separator), or not aggregating this parameter has no effect.

    Returns: outliers
        Returns an OutlierTable object, with outlier and non-outlier counts and metadata
        about how the outliers were called.

    """

    samples = df.columns
    logging.info("Calling outliers for %s samples" % len(samples))

    df = _convert_to_outliers(df, samples, iqrs, up_or_down)
    df = _convert_to_counts(df, samples, aggregate, ind_sep)
    outliers = OutlierTable(df, up_or_down, iqrs, samples, None)

    if save_frac_table:
        frac_path = os.path.abspath(frac_table_file_name % (output_prefix, up_or_down))
        logging.info("Saving outlier fraction table to %s" % frac_path)
        outliers.frac_table.to_csv(frac_path, sep="\t")

    if save_outlier_table:
        out_path = os.path.abspath(
            outlier_table_file_name % (output_prefix, up_or_down)
        )
        logging.info("Saving outlier table to %s" % out_path)
        df.to_csv(out_path, sep="\t")

    return outliers


def compare_groups_outliers(
    outliers: OutlierTable,
    annotations: DataFrame,
    frac_filter: Optional[float] = 0.3,
    save_qvalues: bool = False,
    output_prefix: str = "outliers",
    save_comparison_summaries: bool = False,
) -> qValues:
    """Takes an OutlierTable object and a sample annotation DataFrame and performs comparisons for
    any column in annotations with exactly 2 groups. For each group identified in the annotations
    DataFrame, this function will calculate the q-values of enrichment of outliers for each row in
    each group.

    Args:
        outliers: An OutlierTable, with a DataFrame of outlier and non-outlier counts, \
        as well as parameters for how outliers were calculated.
        annotations: A DataFrame with samples as rows and annotations as columns. Each \
        column must contain exactly 2 different categories, not counting missing values. Columns \
        without 2 options will be ignored.
        frac_filter: The fraction of samples in the group of interest that must \
        have an outlier value to be considered in the comparison. Float between 0 and 1 or None.
        save_qvalues: Whether to write a file with a table of qvalues.
        output_prefix: If files are written, a prefix for the files.
        save_comparison_summaries: Whether to write a file for each annotation column with the \
        counts in the fisher table, pvalues and q values per row.

    Returns: qvals
        A qValues object, which includes a DataFrame of q-values for each comparison, \
        as well as some metadata about how the comparisons were performed.


    """

    df = outliers.df
    samples = outliers.samples
    up_or_down = outliers.up_or_down
    results_df = pd.DataFrame(index=df.index)
    for comp in annotations.columns:
        logging.info("Testing for enrichment in %s comparison" % comp)

        group0_label, group0, group1_label, group1 = get_sample_lists(annotations, comp)
        # Checking everything is in place
        if group0 is None:
            logging.error(
                "There are not exactly 2 groups of samples, skipping %s" % comp
            )
            continue
        not_there = [samp for samp in group0 if samp not in samples] + [
            samp for samp in group1 if samp not in samples
        ]
        if not_there:
            logging.warning(
                "These samples were not found in outliers table: "
                "%s, continuing without them. " % ", ".join(not_there)
            )
        group0 = [samp for samp in group0 if samp in samples]
        group1 = [samp for samp in group1 if samp in samples]
        if len(group0) < 2:
            logging.error(
                "Group %s does not have at least two samples, "
                "skipping comparison %s. " % (group0_label, comp)
            )
            continue
        if len(group1) < 2:
            logging.error(
                "Group %s does not have at least two samples, "
                "skipping comparison%s. " % (group1_label, comp)
            )
            continue

        # doing tests
        label0 = fdr_col_label % (comp, group0_label)
        results_df, fisher_info0 = _compare_groups(
            results_df, df, group0, group1, frac_filter, label0
        )

        label1 = fdr_col_label % (comp, group1_label)
        results_df, fisher_info1 = _compare_groups(
            results_df, df, group1, group0, frac_filter, label1
        )

        if save_comparison_summaries:
            fisher_info0.columns = [
                "%s_%s_%s" % (outlier_count_lab, comp, group0_label),
                "%s_%s_%s" % (outlier_count_lab, comp, group1_label),
                "%s_%s_%s" % (not_outlier_count_lab, comp, group0_label),
                "%s_%s_%s" % (not_outlier_count_lab, comp, group1_label),
                specific_fisher_p % (comp, group0_label),
            ]

            fisher_info1.columns = [
                "%s_%s_%s" % (outlier_count_lab, comp, group1_label),
                "%s_%s_%s" % (outlier_count_lab, comp, group0_label),
                "%s_%s_%s" % (not_outlier_count_lab, comp, group1_label),
                "%s_%s_%s" % (not_outlier_count_lab, comp, group0_label),
                specific_fisher_p % (comp, group1_label),
            ]
            comp_df = pd.concat(
                [fisher_info0, fisher_info1], axis=0, join="outer", sort=True
            ).merge(
                results_df.reindex([label0, label1], axis=1),
                left_index=True,
                right_index=True,
            )
            if len(comp_df) > 0:
                comp_df.to_csv(
                    ind_comparison_file_name % (output_prefix, up_or_down, comp),
                    sep="\t",
                )
    results_df = results_df.dropna(how="all", axis=0)
    if save_qvalues:
        qval_path = os.path.abspath(qvalues_file_name % (output_prefix, up_or_down))
        logging.info("Saving qvalues to %s" % qval_path)
        results_df.to_csv(qval_path, sep="\t")
    qvals = qValues(results_df, annotations.columns, frac_filter)
    return qvals


def deva(
    df: DataFrame,
    annotations: DataFrame,
    iqrs: float = 1.5,
    up_or_down: str = "up",
    aggregate: bool = True,
    save_outlier_table: bool = False,
    save_frac_table: bool = False,
    frac_filter: Optional[float] = 0.3,
    save_qvalues: bool = False,
    output_prefix: str = "outliers",
    ind_sep: str = "-",
    save_comparison_summaries: bool = False,
) -> Tuple[OutlierTable, qValues]:
    """
    Takes a DataFrame of values and returns OutlierTable and qValues objects. This command runs
    the whole outliers pipeline. The DataFrame in the OutlierTable object can be used to run more
    comparisons in future. The qValues object can be used for visualization, or writing
    significant gene lists.

    Args:
        df: Input DataFrame with samples as columns and sites/genes as rows.
        annotations: A DataFrame with samples as rows and annotations as columns. Each \
        column must contain exactly 2 different values, not counting missing \
        values. Other columns will be ignored.
        iqrs: The number of interquartile ranges (IQRs) above or below the median to consider a \
        value as an outlier.
        up_or_down: Whether to call up or down outliers. Up is above the median; down \
        is below the median. Options "up" or "down".
        aggregate: Whether to sum outliers across a grouping (e.g. gene-level) than individual \
        sites. For instance if columns indicate phosphosites on proteins, with the format \
        "RAG2-S365", output will show counts of outliers per protein (e.g. RAG2) rather than on \
        individual sites (e.g. RAG2-S365).
        save_outlier_table: Whether to write a file with the outlier count table.
        save_frac_table: Whether to write a file of the fraction of outliers.
        frac_filter: The fraction of samples in the group of interest that must \
        have an outlier value to be considered in the comparison. Float between 0 and 1 or None.
        save_qvalues: Whether to write a file of qvalues.
        output_prefix: If files are written, a prefix for the files.
        ind_sep: The separator used in the columns, for instance, to separate a gene and site. \
        If just using genes (i.e. no separator), or not aggregating this parameter \
        has no effect.
        save_comparison_summaries: Whether to write a table for each comparison with the \
        counts in the fisher table, pvalues and qvalues per row.

    Returns: outliers, qvals
        Returns an OutlierTable object and qValues object.

    """

    logging.info("Making outliers table")
    outliers = make_outliers_table(
        df,
        iqrs,
        up_or_down,
        aggregate,
        save_outlier_table,
        save_frac_table,
        output_prefix,
        ind_sep,
    )

    logging.info("Performing group comparisons")
    qvals = compare_groups_outliers(
        outliers,
        annotations,
        frac_filter,
        save_qvalues,
        output_prefix,
        save_comparison_summaries,
    )

    return outliers, qvals
