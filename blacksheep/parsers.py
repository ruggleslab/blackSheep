import os.path
import pandas as pd
import numpy as np
from pandas import DataFrame
from typing import Iterable
from blacksheep.classes import OutlierTable
from blacksheep._constants import *


def _is_valid_file(arg: str) -> str:
    """Checks if file exists (probably pretty redundant except as a type checker in argparse)

    Args:
        arg: File path

    Returns: arg
        Validated file path

    """
    if not os.path.exists(arg):
        raise FileNotFoundError("%s does not exist" % arg)
    return arg


def _check_output_prefix(arg: str) -> str:
    """Checks if output prefix is valid

    Args:
        arg: Output prefix

    Returns: arg
        Validation output prefix

    """

    if "/" in arg:
        prefix = arg.rsplit("/", 1)[0]
        _is_valid_file(prefix)
    return arg


def _check_suffix(path: str) -> str:
    """Checks that file is a .csv or .tsv file, and returns which sep to use

    Args:
        path: File path

    Returns: sep
        Sep to use for parsing

    """

    if path[-4:] == ".tsv":
        return "\t"
    if path[-4:] == ".csv":
        return ","
    raise ValueError("File must be .csv or .tsv")


def read_in_values(path: str) -> DataFrame:
    """Figures out sep and parsing file into dataframe.

    Args:
        path: File path

    Returns: df
        DataFrame from table in file

    """
    sep = _check_suffix(path)
    return pd.read_csv(_is_valid_file(path), sep=sep, index_col=0)


def read_in_outliers(path: str, updown: str, iqrs: float) -> OutlierTable:
    """Parses a file into an OutlierTable object.

    Args:
        path: File path
        updown: Whether the outliers represent up or down outliers
        iqrs: How many IQRs were used to define an outlier

    Returns: outliers
        OutlierTable object

    """

    sep = _check_suffix(path)
    df = pd.read_csv(_is_valid_file(path), sep=sep, index_col=0)
    samples = sorted(list(set([ind.rsplit(col_seps, 1)[0] for ind in df.columns])))
    return OutlierTable(df, updown, iqrs, samples, None)


def binarize_annotations(df: DataFrame) -> DataFrame:
    """Takes an annotation DataFrame, checks each column for the number of possible values,
    and adjusts based on that. If the column has 0 or 1 options, it is dropped. Cols with 2
    possible values are retained as-is. Cols with more than 2 values are expanded. For each
    value in that column, a new column is created with val and not_val options.

    Args:
        df: Annotations DataFrame.

    Returns: new_df
        Refactored annotations DataFrame.
    """

    new_df = pd.DataFrame(index=df.index)
    for col in df.columns:
        if len(df[col].dropna().value_counts().keys()) == 2:
            new_df[col] = df[col]
        elif len(df[col].dropna().value_counts().keys()) > 2:
            for val in df[col].dropna().value_counts().keys():
                val = str(val).replace("_", "-")
                new_df.loc[(df[col] != val), binarized_col_name % (col, val)] = (
                    outgroup_val % val
                )
                new_df.loc[(df[col] == val), binarized_col_name % (col, val)] = val
                new_df.loc[df[col].isnull(), binarized_col_name % (col, val)] = np.nan
    return new_df


def normalize(df: DataFrame) -> DataFrame:
    """Performs median of ratios normalization on a given dataframe, then a log2 transform.

    Args:
        df: Unnormalized values dataframe

    Returns: Normalized dataframe

    """
    return np.log2(df.divide(
        df.divide(
            df.mean(axis=1), axis=0).replace([np.inf, -np.inf], np.nan).median()
    ))


def subset_by_genes(
        outliers: DataFrame, ind_list: Iterable[str], ind_sep: str = None,
        ) -> DataFrame:
    if ind_sep:
        return outliers.loc[[i for i in outliers.index if i.split(ind_sep, 1)[0] in ind_list], :]
    return outliers.loc[ind_list, :]