from typing import List, Optional

import pandas as pd
from pandas import DataFrame
from .classes import OutlierTable, qValues
from .outlierTable import convertToOutliers
from .outlierTable import convertToCounts
from .outlierTable import makeFracTable
from .comparisons import compareGroups
from .comparisons import getSampleLists


SampleList = List[str]

def renameFisherTable(string, labelmap):
    for orig, out in labelmap.items():
        string = string.replace(orig, out)
    return string


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
    """

    :param df: Input DataFrame with samples as rows and sites or genes as columns.
    :param iqrs: The number of IQRs above or below the median to consider a value as an outlier.
    Default 1.5.
    :param up_or_down: Whether to call up or down outliers. Up meaning above the median; down
    meaning below the median. Options "up" or "down", default "up".
    :param aggregate: Whether to count up outliers for a more general grouping than individual
    sites. For instance if columns indicate phosphosites on proteins, with the format
    "RAG2-S365", output will show counts of outliers per protein (e.g. RAG2) rather than on
    individual sites (e.g. RAG2-S365). Default True.
    :param frac_table: Whether to output a table representing the fraction of outliers per row,
    per sample. This table is useful for visualization but not downstream analysis. Default False.
    :param ind_sep: The separator used in the columns, for instance, to separate a gene and site.
    If just using genes (i.e. no separator), or not aggregating this parameter has no effect.
    Default "-"
    :param save_outlier_table: Whether to write a file with the outlier count table. Default False.
    :param save_frac_table: Whether to write a file with the outlier fraction table. Default False.
    :param output_prefix: If files are written, a prefix for the files.
    :return: Always returns a outliers count table, with outlier and non-outlier counts of
    samples as rows and genes/sites as columns. If frac_table is set to True, also returns a
    table of fractions of outliers with samples as rows, and genes/sites as columns.
    """
    df = df.transpose()
    samples = df.columns

    df = convertToOutliers(df, samples, iqrs, up_or_down)

    df = convertToCounts(df, samples, aggregate, ind_sep).transpose()

    fracTable = makeFracTable(df, samples).transpose().dropna(how='all')
    if save_frac_table:
        fracTable.to_csv(
            "%s.%s.fraction_table.tsv" % (output_prefix, up_or_down), sep="\t"
        )

    if save_outlier_table:
        df.to_csv(
            "%s.%s.count_table.tsv" % (output_prefix, up_or_down), sep="\t"
        )

    outliers = OutlierTable(df, up_or_down, iqrs, samples, fracTable)
    return outliers


def compare_groups_outliers(
    outliers: OutlierTable,
    annotations: DataFrame,
    frac_filter: Optional[float] = 0.3,
    save_qvalues: bool = False,
    output_prefix: str = "outliers",
    output_comparison_summaries: bool = False,
) -> qValues:
    """

    :param outliers:
    :param annotations: A DataFrame with samples as the index and annotations as columns. Each
    column must contain exactly 2 different values, and optionally missing values. Columns with
    less or more than 2 options will be ignored.
    :param frac_filter: The fraction of samples in group0 (i.e. the group of interest) that must
    have an outlier value to be considered in the comparison. Float between 0 and 1 or None.
    :param save_qvalues: Whether to output a table of qvalues. Default False.
    :param output_prefix: If files are written, a prefix for the files. Default "outliers"
    :param up_or_down: Whether the input file is up or down outliers. Default "up"
    :param output_comparison_summaries: Whether to write a table for each comparison with the
    counts in the fisher table, pvalues and q values per row. Default False.
    :return:
    """
    # TODO comvert prints to logging
    df = outliers.df.transpose()
    outlier_samples = outliers.samples
    up_or_down = outliers.up_or_down
    results_df = pd.DataFrame(index=df.index)
    for comp in annotations.columns:

        group0_label, group0, group1_label, group1 = getSampleLists(annotations, comp)
        # Checking everything is in place
        if group0 is None:
            print("Number of categories in %s is not 2, skipping %s" % (comp, comp))
            continue
        not_there = [samp for samp in group0 if not samp in outlier_samples] + [
            samp for samp in group1 if not samp in outlier_samples
        ]
        group0 = [samp for samp in group0 if samp in outlier_samples]
        group1 = [samp for samp in group1 if samp in outlier_samples]
        if (len(group0) < 2) or (len(group1) < 2):
            print(
                "Number of categories in %s is not 2, or too few samples in each category, "
                "skipping %s" % (comp, comp)
            )
            continue
        if len(not_there) > 0:
            print(
                "Samples %s missing in outlier table, continuing %s without them"
                % (", ".join(not_there), comp)
            )

        # doing tests
        label0 = "fisherFDR_%s_%s" % (comp, group0_label)
        results_df, fisher_info0 = compareGroups(
            results_df, df, group0, group1, frac_filter, label0
        )

        label1 = "fisherFDR_%s_%s" % (comp, group1_label)
        results_df, fisher_info1 = compareGroups(
            results_df, df, group1, group0, frac_filter, label1
        )

        if output_comparison_summaries:
            label_map = {
                "0": "_%s_%s" % (comp, group0_label),
                "1": "_%s_%s" % (comp, group1_label),
                "fisherp": "enrichment_fisherp_%s_%s" % (comp, group0_label),
            }
            fisher_info0.columns = [
                renameFisherTable(col, label_map) for col in \
                    fisher_info0.columns
            ]

            label_map = {
                "0": "_%s_%s" % (comp, group1_label),
                "1": "_%s_%s" % (comp, group0_label),
                "fisherp": "enrichment_fisherp_%s_%s" % (comp, group1_label),
            }
            fisher_info1.columns = [
                renameFisherTable(col, label_map) for col in fisher_info1.columns
            ]

            comp_df = pd.concat(
                [fisher_info0, fisher_info1], axis=0,join="outer",
                sort=True
            ).merge(results_df[
                [label0, label1]
            ], left_index=True, right_index=True)
            comp_df.to_csv(
                "%s.%s.%s.qvalues.tsv" % (output_prefix, up_or_down, comp), sep="\t"
            )

    if save_qvalues:
        results_df.to_csv("%s.%s.qvalues.tsv" % (output_prefix, up_or_down), sep="\t")
    qvals = qValues(df, annotations.columns, frac_filter)
    return qvals


def run_outliers(
    df: DataFrame,
    annotations: DataFrame,
    frac_filter: Optional[float] = 0.3,
    iqrs: float = 1.5,
    up_or_down: str = "up",
    aggregate: bool = True,
    save_outlier_table: bool = False,
    save_frac_table: bool = False,
    save_qvalues: bool = False,
    output_prefix: str = "outliers",
    ind_sep: str = "-",
    output_comparison_summaries: bool = False,
) -> DataFrame:
    """

    :param df: Input DataFrame with samples as rows and sites or genes as columns.
    :param annotations: A DataFrame with samples as the index and annotations as columns. Each
    column must contain exactly 2 different values, and optionally missing values. Columns with
    less or more than 2 options will be ignored.
    :param frac_filter: The fraction of samples in group0 (i.e. the group of interest) that must
    have an outlier value to be considered in the comparison. Float between 0 and 1 or None.
    :param iqrs: The number of IQRs above or below the median to consider a value as an outlier.
    Default 1.5.
    :param up_or_down: Whether to call up or down outliers. Up meaning above the median; down
    meaning below the median. Options "up" or "down", default "up".
    :param aggregate: Whether to count up outliers for a more general grouping than individual
    sites. For instance if columns indicate phosphosites on proteins, with the format
    "RAG2-S365", output will show counts of outliers per protein (e.g. RAG2) rather than on
    individual sites (e.g. RAG2-S365). Default True.
    :param frac_table:
    :param save_outlier_table: Whether to write a file with the outlier count table. Default False.
    per row, per sample. This table is useful for visualization but not downstream analysis.
    Default False.
    :param save_frac_table: Whether to output a table representing the fraction of outliers.
    Default False.
    :param save_qvalues: Whether to output a table of qvalues. Default False.
    :param output_prefix: If files are written, a prefix for the files. Default "outliers"
    :param ind_sep: The separator used in the columns, for instance, to separate a gene and site.
    If just using genes (i.e. no separator), or not aggregating this parameter has no effect.
    Default "-"
    :param output_comparison_summaries: Whether to write a table for each comparison with the
    counts in the fisher table, pvalues and q values per row. Default False.
    :return: Always returns a outliers count table, with outlier and non-outlier counts of
    samples as rows and genes/sites as columns. If frac_table is set to True, also returns a
    table of fractions of outliers with samples as rows, and genes/sites as columns.
    """
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

    frac_table = outliers.frac_table
    qvals = compare_groups_outliers(
        outliers,
        annotations,
        frac_filter,
        save_qvalues,
        output_prefix,
        output_comparison_summaries,
    )

    return qvals, frac_table
