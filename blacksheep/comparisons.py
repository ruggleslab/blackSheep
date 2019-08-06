import logging
from typing import List, Tuple, Iterable, Optional
import pandas as pd
from pandas import DataFrame
from pandas import Series
import numpy as np
import scipy.stats
from blacksheep.constants import *


SampleList = List[str]
logger = logging.getLogger("cli")


def _multi_hyp_correct(
    pvalues: Iterable[float], correction_type: str = "Benjamini-Hochberg"
) -> Iterable[float]:
    """Corrects p-values for multiple hypothesis testing

    Args:
        pvalues:  Array of p-values to correct
        correction_type: correction_type: Which procedure to use. Options are "Benjamini-Hochberg" or
        "Bonferroni"

    Returns: rray of p-values corrected for multiple hypothesis testing (aka q-values).

    """

    pvalues = np.array(pvalues)
    n = sum(~np.isnan(pvalues))
    new_pvalues = np.empty(len(pvalues))

    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues

    elif correction_type == "Benjamini-Hochberg":

        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values = [x for x in values if ~np.isnan(x[0])]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n / rank) * pvalue)
        for i in range(0, int(n) - 1):
            if new_values[i] < new_values[i + 1]:
                new_values[i + 1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]

    new_pvalues[np.isnan(pvalues)] = np.nan
    for i, val in enumerate(new_pvalues):
        if ~np.isnan(val):
            new_pvalues[i] = min(1, val)

    return new_pvalues


def get_sample_lists(
    annotations: DataFrame, col: str
) -> Tuple[Optional[str], Optional[SampleList], Optional[str], Optional[SampleList]]:
    """Finds groupings of samples from an annotation DataFrame column.

    Args:
        annotations: A DataFrame with samples as the index and annotations as columns. Each
        column must contain exactly 2 different values, and optionally missing values. Columns with
        less or more than 2 options will be ignored.
        col: Which column for which to define groups.

    Returns: A label for group0, the list of samples in group0, a label for group1 and the list
        of samples in group1.

    """

    groups = list(pd.Series(annotations[col].value_counts().keys()).dropna())
    if len(groups) != 2:
        return None, None, None, None
    group0 = list(annotations.loc[annotations[col] == groups[0], :].index)
    group1 = list(annotations.loc[annotations[col] == groups[1], :].index)
    return groups[0], group0, groups[1], group1


def _filter_outliers(
    df: DataFrame,
    group0_list: SampleList,
    group1_list: SampleList,
    frac_filter: Optional[float],
) -> DataFrame:
    """Filters an outlier count table for rows that are enriched for outliers in group0 and that
    have more than a frac_filter fraction of samples of group0 with an outlier.

    Args:
        df: Outliers count table, output from convertToCounts. Samples are columns,
        genes/sites are the index.
        group0_list: List of samples in the group of interest.
        group1_list: List of samples in the outgroup.
        frac_filter: The fraction of samples in group0 (i.e. the group of interest) that must
        have an outlier value to be considered in the comparison. Float between 0 and 1 or None.

    Returns:A DataFrame with rows that are not enriched in group0 removed. If frac_filter > 0,
    rows without enough outliers in group0 are also removed.

    """

    group0_outliers = [x + col_seps + col_outlier_suffix for x in group0_list]
    group0_notOutliers = [x + col_seps + col_not_outlier_suffix for x in group0_list]
    group1_outliers = [x + col_seps + col_outlier_suffix for x in group1_list]
    group1_notOutliers = [x + col_seps + col_not_outlier_suffix for x in group1_list]

    if frac_filter != None:
        min_num_outlier_samps = len(group0_list) * frac_filter
        num_outlier_samps = (df[group0_outliers] > 0).sum(axis=1)
        df = df.loc[num_outlier_samps >= min_num_outlier_samps, :]

    # Filter for higher proportion of outliers in group0 than group1
    group0_outlier_rate = (
        df[group0_outliers]
        .sum(axis=1)
        .divide(df[group0_outliers + group0_notOutliers].sum(axis=1), axis=0)
    )
    group1_outlier_rate = (
        df[group1_outliers]
        .sum(axis=1)
        .divide(df[group0_outliers + group1_notOutliers].sum(axis=1), axis=0)
    )

    df = df.loc[group0_outlier_rate > group1_outlier_rate, :]

    return df


def _fisher_test_groups(
    group0_list: SampleList,
    group1_list: SampleList,
    outlier_table: DataFrame,
    correction_type: str = mult_hypoth_method,
) -> Tuple[Series, DataFrame]:
    """Performs fishers test by counting outlier and not outlier sites in two groups. Corrects for
    multiple hypothesis testing.

    Args:
        group0_list: List of samples in group of interest
        group1_list: List of samples in outgroup
        outlier_table: Outlier count table, like output of convertToCounts
        correction_type: Method to use for multiple hypothesis correction.

    Returns: Series of qvalues with index matching filtered rows.

    """

    outliers_group0_list = [x + col_seps + col_outlier_suffix for x in group0_list]
    notOutliers_group0_list = [
        x + col_seps + col_not_outlier_suffix for x in group0_list
    ]
    outliers_group1_list = [x + col_seps + col_outlier_suffix for x in group1_list]
    notOutliers_group1_list = [
        x + col_seps + col_not_outlier_suffix for x in group1_list
    ]

    outlier_table[outlier_count_lab + general_group_label_0] = outlier_table[
        outliers_group0_list
    ].sum(axis=1)
    outlier_table[not_outlier_count_lab + general_group_label_0] = outlier_table[
        notOutliers_group0_list
    ].sum(axis=1)
    outlier_table[outlier_count_lab + general_group_label_1] = outlier_table[
        outliers_group1_list
    ].sum(axis=1)
    outlier_table[not_outlier_count_lab + general_group_label_1] = outlier_table[
        notOutliers_group1_list
    ].sum(axis=1)

    outlier_table[fisherp_col] = outlier_table.apply(
        (
            lambda r: scipy.stats.fisher_exact(
                [
                    [
                        r[outlier_count_lab + general_group_label_0],
                        r[outlier_count_lab + general_group_label_1],
                    ],
                    [
                        r[not_outlier_count_lab + general_group_label_0],
                        r[not_outlier_count_lab + general_group_label_1],
                    ],
                ]
            )[1]
        ),
        axis=1,
    )

    outlier_table[fisherfdr_col] = _multi_hyp_correct(
        list(outlier_table[fisherp_col]), correction_type
    )
    fisher_info = outlier_table[
        [
            outlier_count_lab + general_group_label_0,
            outlier_count_lab + general_group_label_1,
            not_outlier_count_lab + general_group_label_0,
            not_outlier_count_lab + general_group_label_1,
            fisherp_col,
        ]
    ]
    return outlier_table[fisherfdr_col], fisher_info


def _compare_groups(
    results_df: DataFrame,
    outliers: DataFrame,
    group0: SampleList,
    group1: SampleList,
    frac_filter: Optional[float],
    label: str,
) -> Tuple[DataFrame, Optional[DataFrame]]:
    """Performs fisher test and cleans up a fisher infor table for making output for each comparison

    Args:
        results_df: Accumulating qvalues DataFrame
        outliers: Outliers count DataFrame
        group0: List of samples in group of interest
        group1: List of samples in outgroup
        frac_filter: Fraction of samples in group of interest require to have an outlier per
    site to be considered in analysis
        label: What to call the FDR output column on the qvalues DataFrame

    Returns: Concatenated qvalues DataFrame and a table of info about the comparison

    """

    df = _filter_outliers(outliers, group0, group1, frac_filter)
    logger.info("Calculating enrichment in %s rows for %s" % (len(df), label))
    if len(df) > 0:
        col, fisher_info = _fisher_test_groups(group0, group1, df)
        col = DataFrame(col)
        col.columns = [label]
        results_df = pd.concat([results_df, col], axis=1, join="outer", sort=False)
    else:
        logger.warning("No rows tested for %s" % label)
        fisher_info = DataFrame(columns=[
            outlier_count_lab + general_group_label_0,
            outlier_count_lab + general_group_label_1,
            not_outlier_count_lab + general_group_label_0,
            not_outlier_count_lab + general_group_label_1,
            fisherp_col,
        ])
    return results_df, fisher_info
