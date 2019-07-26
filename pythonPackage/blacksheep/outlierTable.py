import pandas as pd
import numpy as np
import scipy.stats
from pandas import DataFrame
from typing import List

SampleList = List[str]


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
    df = df.transpose()
    cols_outliers = [x + "_outliers" for x in samples]
    cols_notOutliers = [x + "_notOutliers" for x in samples]
    df = df.fillna(0)
    num_total_psites = df[cols_notOutliers].values + df[cols_outliers].values
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_outliers = df[cols_outliers].values / num_total_psites

    frac_outliers = pd.DataFrame(
        frac_outliers, index=df.index, columns=samples
    ).transpose()

    return frac_outliers
