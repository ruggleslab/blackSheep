import pandas as pd
import numpy as np
import scipy.stats
from pandas import DataFrame
from pandas import Series
from typing import List, Tuple, Iterable, Optional

SampleList = List[str]


def multHypothCorrect(
    pvalues: Iterable[float], correction_type: str = "Benjamini-Hochberg"
) -> Iterable[float]:
    """

    :param pvalues: Array of p-values to correct
    :param correction_type: Which procedure to use. Options are "Benjamini-Hochberg" or
    "Bonferroni"; Default is "Benjamini-Hochberg"
    :return: Array of p-values corrected for multiple hypothesis testing (aka q-values).
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


def getSampleLists(
    annotations: DataFrame, col: str
) -> Tuple[Optional[str], Optional[SampleList], Optional[str], Optional[SampleList]]:
    """

    :param annotations: A DataFrame with samples as the index and annotations as columns. Each
    column must contain exactly 2 different values, and optionally missing values. Columns with
    less or more than 2 options will be ignored.
    :param col: Which column to define groups for.
    :return: A label for group0, the list of samples in group0, a label for group1 and the list
    of samples in group1.
    """
    groups = list(pd.Series(annotations[col].value_counts().keys()).dropna())
    if len(groups) != 2:
        return None, None, None, None
    group0 = list(annotations.loc[annotations[col] == groups[0], :].index)
    group1 = list(annotations.loc[annotations[col] == groups[1], :].index)
    return groups[0], group0, groups[1], group1


def filterOutliers(
    df: DataFrame,
    group0_list: SampleList,
    group1_list: SampleList,
    frac_filter: Optional[float],
) -> DataFrame:
    """

    :param df: Outliers count table, output from convertToCounts. Samples are columns,
    genes/sites are the index.
    :param group0_list: List of samples in the group of interest.
    :param group1_list: List of samples in the outgroup.
    :param frac_filter: The fraction of samples in group0 (i.e. the group of interest) that must
    have an outlier value to be considered in the comparison. Float between 0 and 1 or None.
    :return: A DataFrame with rows that are not enriched in group0 removed. If frac_filter > 0,
    rows without enough outliers in group0 are also removed.
    """

    group0_outliers = [x + "_outliers" for x in group0_list]
    group0_notOutliers = [x + "_notOutliers" for x in group0_list]
    group1_outliers = [x + "_outliers" for x in group1_list]
    group1_notOutliers = [x + "_notOutliers" for x in group1_list]

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


def testDifferentGroupsOutliers(
    group0_list: SampleList,
    group1_list: SampleList,
    outlier_table: DataFrame,
    correction_type: str = "Benjamini-Hochberg",
) -> Series:
    """

    :param group0_list: List of samples in group of interest
    :param group1_list: List of samples in outgroup
    :param outlier_table: Outlier count table, like output of convertToCounts
    :return: Series of qvalues with index matching filtered rows.
    """

    outliers_group0_list = [x + "_outliers" for x in group0_list]
    notOutliers_group0_list = [x + "_notOutliers" for x in group0_list]
    outliers_group1_list = [x + "_outliers" for x in group1_list]
    notOutliers_group1_list = [x + "_notOutliers" for x in group1_list]

    outlier_table["Outlier0"] = outlier_table[outliers_group0_list].sum(axis=1)
    outlier_table["NotOutlier0"] = outlier_table[notOutliers_group0_list].sum(axis=1)
    outlier_table["Outlier1"] = outlier_table[outliers_group1_list].sum(axis=1)
    outlier_table["NotOutlier1"] = outlier_table[notOutliers_group1_list].sum(axis=1)

    outlier_table["fisherp"] = outlier_table.apply(
        (
            lambda r: scipy.stats.fisher_exact(
                [[r["Outlier0"], r["Outlier1"]], [r["NotOutlier0"], r["NotOutlier1"]]]
            )[1]
        ),
        axis=1,
    )

    outlier_table["fisherFDR"] = multHypothCorrect(
        list(outlier_table["fisherp"]), correction_type
    )
    fisher_info = outlier_table[
        ["Outlier0", "NotOutlier0", "Outlier1", "NotOutlier1", "fisherp"]
    ]
    return outlier_table["fisherFDR"], fisher_info


def compareGroups(
    results_df: DataFrame,
    outliers: DataFrame,
    group0: SampleList,
    group1: SampleList,
    frac_filter: Optional[float],
    label: str,
) -> Tuple[DataFrame, Optional[DataFrame]]:
    df = filterOutliers(outliers, group0, group1, frac_filter)
    if len(df) > 0:
        print(
            "Testing %s rows for enrichment in %s %s samples"
            % (len(df), comp, group0_label)
        )
        col, fisher_info = testDifferentGroupsOutliers(group0, group1, df)
        col = pd.DataFrame(col)
        col.columns = [label]
        results_df = pd.concat([results_df, col], axis=1, join="outer", sort=False)
    else:
        print(
            "No rows had outliers in at least %s of %s %s samples"
            % (frac_filter, comp, group0_label)
        )
    return results_df, fisher_info
