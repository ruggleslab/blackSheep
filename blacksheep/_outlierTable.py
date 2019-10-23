import pandas as pd
import numpy as np
import scipy.stats
from pandas import DataFrame
from typing import List
from blacksheep._constants import *


SampleList = List[str]


def _convert_to_outliers(
    df: DataFrame, samples: SampleList, num_iqrs: float, up_or_down: str
) -> DataFrame:
    """Calls outliers on a given values table.

    Args:
        df: Input DataFrame with samples as columns and genes or sites as rows. \
        Index should be an identifier for each row.
        samples: List of samples to be considered in the distribution when defining median, \
        IQR and outliers.
        num_iqrs: How many inter-quartile ranges (IQRs) above or below the median to consider \
        something an outlier.
        up_or_down: Whether to call outliers above the median (up) or below the median (down)

    Returns:
        A DataFrame with outlier calls for each value. 0 means not an outlier; 1 means there is \
        an outlier. Missing values are propagated.

    """
    if num_iqrs <= 0:
        raise ValueError("num_iqrs must be greater than 0")

    df = df.copy()
    df[row_iqr_name] = scipy.stats.iqr(df[samples], axis=1, nan_policy="omit")
    df[row_median_name] = np.nanquantile(df[samples], q=0.5, axis=1)

    outlier_df = pd.DataFrame()

    if up_or_down == "up":
        df[row_upper_bound_name] = df[row_median_name] + (num_iqrs * df[row_iqr_name])
        outlier_df[samples] = (
            df[samples].gt(df[row_upper_bound_name], axis=0).astype(int)
        )
        outlier_df[df[samples].isnull()] = np.nan
        return outlier_df
    elif up_or_down == "down":
        df[row_lower_bound_name] = df[row_median_name] - (num_iqrs * df[row_iqr_name])
        outlier_df[samples] = (
            df[samples].lt(df[row_lower_bound_name], axis=0).astype(int)
        )
        outlier_df[df[samples].isnull()] = np.nan
        return outlier_df
    raise ValueError("up_or_down must be either 'up' or 'down'")


def _convert_to_counts(
    df: DataFrame, samples: SampleList, aggregate: bool, ind_sep: str
) -> DataFrame:
    """Counts outliers and non-outlier values for each sample and each row (if aggregate=False)
    or each unique identifier (if aggregate=True).

    Args:
        df: Outlier DataFrame from convertToOutliers function.
        samples: List of samples to consider. Should be same list as input to convertToOutliers.
        aggregate: Whether or not to collapse multiple rows with identical identifiers before \
        a separater. Collapsing values is done by counting outliers and non-outliers for each \
        sample per unique identifier.
        ind_sep: The separator used in the index to separate a more general ID and less specific \
        ID. e.g. in RAG2-S365 the separater is "-".

    Returns:
        Outlier DataFrame that has a column with counts of non-outliers and outliers for \
        each sample, with a row for each input site (if no aggregation) or each unique identifier \
        (with aggregation). This table is the input for the comparison function.

    """
    not_outlier_cols = [x + col_seps + col_not_outlier_suffix for x in samples]
    outlier_cols = [x + col_seps + col_outlier_suffix for x in samples]

    if aggregate:
        df.index = [ind.split(ind_sep)[0] for ind in df.index]
        df_inv = df == 0

        output_df = pd.DataFrame()
        output_df[not_outlier_cols] = df_inv.groupby(level=0)[samples].sum()
        output_df[outlier_cols] = df.groupby(level=0)[samples].sum()
    elif not aggregate:
        output_df = pd.DataFrame(index=df.index)
        output_df[outlier_cols] = df[samples]
        output_df[not_outlier_cols] = 1 - df[samples]
    return output_df
