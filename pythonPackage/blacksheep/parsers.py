import os.path
from typing import Iterable
import pandas as pd
from pandas import DataFrame
import numpy as np
from .classes import OutlierTable
from .constants import *


def is_valid_file(arg: str) -> str:
    """
    Checks if file exists (probably pretty redundant except as a type checker in argparse)
    :param arg: File path
    :return: Validated file path
    """
    if not os.path.exists(arg):
        raise FileNotFoundError("%s does not exist" % arg)
    return arg


def check_output_prefix(arg: str) -> str:
    """
    Checks if output prefix is valid
    :param arg: Output prefix
    :return: Validation output prefix
    """
    if "/" in arg:
        prefix = arg.rsplit("/", 1)[0]
        is_valid_file(prefix)
    return arg


def check_suffix(path: str) -> str:
    """
    Checks that file is a .csv or .tsv file, and returns which sep to use
    :param path: File path
    :return: Sep to use for parsing
    """
    if path[-4:] == ".tsv":
        sep = "\t"
    elif path[-4:] == ".csv":
        sep = ","
    else:
        raise ValueError("File must be .csv or .tsv")
    return sep


def parse_values(path: str) -> DataFrame:
    """
    Figures out sep and parsing file into dataframe.
    :param path: File path
    :return: DataFrame from table in file
    """
    sep = check_suffix(path)
    df = pd.read_csv(is_valid_file(path), sep=sep, index_col=0)
    return df


def parse_outliers(path: str, updown: str, iqrs: float) -> OutlierTable:
    """
    Parses a file into an OutlierTable object.
    :param path: File path
    :param updown: Whether the outliers represent up or down outliers
    :param iqrs: How many IQRs were used to define an outlier
    :return: OutlierTable object
    """
    sep = check_suffix(path)
    df = pd.read_csv(is_valid_file(path), sep=sep, index_col=0)
    samples = sorted(list(set([ind.rsplit(col_seps, 1)[0] for ind in df.columns])))
    outliers = OutlierTable(df, updown, iqrs, samples, None)
    outliers.make_frac_table()
    return outliers


def binarize_annotations(df: DataFrame) -> DataFrame:
    """
    Takes an annotation DataFrame, checks each column for the number of possible values,
    and adjusts based on that. If the column has 0 or 1 options, it is dropped. Cols with 2
    possible values are retained as-is. Cols with more than 2 values are expanded. For each
    value in that column, a new column is created with val and not_val options.
    :param df: Annotations DataFrame.
    :return: Refactored annotations DataFrame.
    """
    new_df = pd.DataFrame(index=df.index)
    for col in df.columns:
        if len(df[col].dropna().value_counts().keys()) == 2:
            new_df[col] = df[col]
        elif len(df[col].dropna().value_counts().keys()) > 2:
            for val in df[col].dropna().value_counts().keys():
                new_df.loc[(df[col] != val), binarized_col_name % (col, val)] = (
                    outgroup_val % val
                )
                new_df.loc[(df[col] == val), binarized_col_name % (col, val)] = val
                new_df.loc[df[col].isnull(), binarized_col_name % (col, val)] = np.nan
    return new_df


def list_to_file(lis: Iterable, filename: str):
    """
    Takes an iterable and a file path and writes a value per line from the iterable into the new
    file.
    :param lis: Iterable to write to file
    :param filename: Filename to write to.
    :return: None
    """
    with open(check_output_prefix(filename), "w") as fh:
        for x in lis:
            fh.write("%s\n" % x)
