import pandas as pd
import numpy as np
import scipy.stats
from pandas import DataFrame
from pandas import Series
from typing import List, Tuple, Iterable, Optional

SampleList = List[str]

# Converting to outliers
def convertToOutliers(
    df: DataFrame, samples: SampleList, NUM_IQRs: float, up_or_down: str
) -> DataFrame:
    """
    :param df: Input DataFrame with samples as columns and genes or sites as rows.
    Index should be an identifier for each row.
    :param samples: List of samples to be considered in the distribution when defining median,
    IQR and outliers.
    :param NUM_IQRs: How many inter-quartile ranges (IQRs) above or below the median to consider
    something an outlier.
    :param up_or_down: Whether to call outliers above the median (up) or below the median (down)
    :return: A DataFrame with outlier calls for each value. 0 means not an outlier; 1 means there
    is an outlier.Missing values are propagated.
    """
    # TODO add threshold outputs
    df["row_iqr"] = scipy.stats.iqr(df[samples], axis=1, nan_policy="omit")
    df["row_median"] = np.nanquantile(df[samples], q=0.5, axis=1)

    outlier_df = pd.DataFrame()

    if up_or_down == "up":
        df["row_medPlus"] = df["row_median"] + (NUM_IQRs * df["row_iqr"])
        outlier_df[samples] = df[samples].gt(df["row_medPlus"], axis=0).astype(int)

    elif up_or_down == "down":
        df["row_medMinus"] = df["row_median"] - (NUM_IQRs * df["row_iqr"])
        outlier_df[samples] = df[samples].lt(df["row_medMinus"], axis=0).astype(int)

    outlier_df[df[samples].isnull()] = np.nan

    return outlier_df


def convertToCounts(
    df: DataFrame, samples: SampleList, aggregate: bool, ind_sep: str
) -> DataFrame:
    """

    :param df: Outlier DataFrame from convertToOutliers function.
    :param samples: List of samples to consider. Should be same list as input to convertToOutliers.
    :param aggregate: Whether or not to collapse multiple rows with identical identifiers before
    a separater. Collapsing values is done by counting outliers and non-outliers for each sample
    per unique identifier.
    :param ind_sep: The separator used in the index to separate a more general ID and less specific
    ID. e.g. in RAG2-S365 the separater is "-".
    :return: Outlier DataFrame that has a column with counts of non-outliers and outliers for
    each sample, with a row for each input site (if no aggregation) or each unique identifier (
    with aggregation).
    This table is the input for the comparison function.
    """
    not_outlier_cols = [x + "_notOutliers" for x in samples]
    outlier_cols = [x + "_outliers" for x in samples]

    if aggregate:
        agg_col = "gene"
        df[agg_col] = [ind.split(ind_sep)[0] for ind in df.index]
        output_df = pd.DataFrame()
        output_df[not_outlier_cols] = df.groupby(by=agg_col)[samples].agg(
            lambda x: pd.Series(x == 0).sum()
        )
        output_df[outlier_cols] = df.groupby(by=agg_col)[samples].agg(
            lambda x: pd.Series(x == 1).sum()
        )
    elif not aggregate:
        output_df = pd.DataFrame(index=df.index)
        output_df[outlier_cols] = df[samples]
        output_df[not_outlier_cols] = 1 - df[samples]
    return output_df


def makeFracTable(df: DataFrame, samples: SampleList):
    """

    :param df: The outlier table that is output from convertToCounts
    :param samples: List of samples input to convertToOutliers and convertToCounts
    :return: A table with one column per sample with the fraction of sites in each row
    that are outliers per each sample. This table is useful for visualization but not statistics.
    """
    cols_outliers = [x + "_outliers" for x in samples]
    cols_notOutliers = [x + "_notOutliers" for x in samples]
    df = df.fillna(0)
    num_total_psites = df[cols_notOutliers].values + df[cols_outliers].values
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_outliers = df[cols_outliers].values / num_total_psites

    frac_outliers = pd.DataFrame(frac_outliers, index=df.index, columns=samples)

    return frac_outliers


def make_outliers_table(
    df: DataFrame,
    iqrs: float = 1.5,
    up_or_down: str = "up",
    aggregate: bool = True,
    frac_table: bool = False,
    save_outlier_table: bool = False,
    save_frac_table: bool = False,
    output_prefix: str = "outliers",
    ind_sep: str = "-",
) -> Tuple[DataFrame, Optional[DataFrame]]:
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
    # TODO add a option for outputing thresholds for outliers
    df = df.transpose()
    samples = df.columns

    df = convertToOutliers(df, samples, iqrs, up_or_down)
    df = convertToCounts(df, samples, aggregate, ind_sep)

    if frac_table:
        fracTable = makeFracTable(df, samples)
        if save_frac_table:
            fracTable.transpose().to_csv(
                "%s.fraction_table.tsv" % output_prefix, sep="\t"
            )
        return df.transpose(), fracTable.transpose()
    if save_outlier_table:
        df.transpose().to_csv("%s.count_table.tsv" % output_prefix, sep="\t")
    return df.transpose(), None


# Making comparisons
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
    df: DataFrame, group0_list: SampleList, group1_list: SampleList, frac_filter: Optional[float]
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
    return outlier_table["fisherFDR"]


def compareGroups(
    results_df: DataFrame,
    outliers: DataFrame,
    group0: SampleList,
    group1: SampleList,
    frac_filter: Optional[float],
    label: str,
) -> DataFrame:
    df = filterOutliers(outliers, group0, group1, frac_filter)
    if len(df) > 0:
        print(
            "Testing %s rows for enrichment in %s %s samples"
            % (len(df), comp, group0_label)
        )
        col = pd.DataFrame(testDifferentGroupsOutliers(group0, group1, df))
        col.columns = [label]
        results_df = pd.concat([results_df, col], axis=1, join="outer", sort=False)
    else:
        print(
            "No rows had outliers in at least %s of %s %s samples"
            % (frac_filter, comp, group0_label)
        )
    return results_df


def compare_groups_outliers(
    outliers: DataFrame,
    annotations: DataFrame,
    frac_filter: Optional[float] = 0.3,
    save_qvalues: bool = False,
    output_prefix: str = "outliers",
) -> DataFrame:
    """

    :param outliers:
    :param annotations: A DataFrame with samples as the index and annotations as columns. Each
    column must contain exactly 2 different values, and optionally missing values. Columns with
    less or more than 2 options will be ignored.
    :param frac_filter: The fraction of samples in group0 (i.e. the group of interest) that must
    have an outlier value to be considered in the comparison. Float between 0 and 1 or None.
    :param save_qvalues:
    :param output_prefix:
    :return:
    """
    outliers = outliers.transpose()
    results_df = pd.DataFrame(index=outliers.index)
    for comp in annotations.columns:
        group0_label, group0, group1_label, group1 = getSampleLists(annotations, comp)
        if group0_label is None:
            print("Number of categories in %s is not 2, skipping %s" % (comp, comp))
            continue

        outlier_samples = [col.rsplit("_", 1)[0] for col in outliers.columns]
        not_there = [samp for samp in group0 if not samp in outlier_samples] + [
            samp for samp in group1 if not samp in outlier_samples
        ]
        if len(not_there) > 0:
            print(
                "Samples %s missing in outlier table, continuing %s without them"
                % (", ".join(not_there), comp)
            )

        group0 = [samp for samp in group0 if samp in outlier_samples]
        group1 = [samp for samp in group1 if samp in outlier_samples]

        if (len(group0) > 0) & (len(group1) > 0):
            label = "%s_%s_enrichment_FDR" % (comp, group0_label)
            results_df = compareGroups(
                results_df, outliers, group0, group1, frac_filter, label
            )

            label = "%s_%s_enrichment_FDR" % (comp, group1_label)
            results_df = compareGroups(
                results_df, outliers, group1, group0, frac_filter, label
            )
        elif len(group0) > 0:
            print("No samples from %s %s in outlier table" % (comp, group0_label))
        elif len(group1) > 0:
            print("No samples from %s %s in outlier table" % (comp, group1_label))
    if save_qvalues:
        results_df.to_csv("%s.qvalues.tsv" % (output_prefix), sep="\t")
    return results_df


def run_outliers(
    df: DataFrame,
    annotations: DataFrame,
    frac_filter: Optional[float] = 0.3,
    iqrs: float = 1.5,
    up_or_down: str = "up",
    aggregate: bool = True,
    frac_table: bool = False,
    save_outlier_table: bool = False,
    save_frac_table: bool = False,
    save_qvalues: bool = False,
    output_prefix: str = "outliers",
    ind_sep: str = "-",
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
    :return: Always returns a outliers count table, with outlier and non-outlier counts of
    samples as rows and genes/sites as columns. If frac_table is set to True, also returns a
    table of fractions of outliers with samples as rows, and genes/sites as columns.
    """
    counts, frac_table = make_outliers_table(
        df,
        iqrs,
        up_or_down,
        aggregate,
        frac_table,
        save_outlier_table,
        save_frac_table,
        output_prefix,
        ind_sep,
    )
    qvalues = compare_groups_outliers(
        counts, annotations, frac_filter, save_qvalues, output_prefix
    )

    return qvalues, frac_table
