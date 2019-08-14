import warnings
import os.path
from typing import Iterable, Optional
import pandas as pd
from pandas import DataFrame, Series
import numpy as np
import sklearn.linear_model as lm
from sklearn.exceptions import UndefinedMetricWarning

from deva.classes import OutlierTable
from deva.constants import *


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
        sep = "\t"
    elif path[-4:] == ".csv":
        sep = ","
    else:
        raise ValueError("File must be .csv or .tsv")
    return sep


def parse_values(path: str) -> DataFrame:
    """Figures out sep and parsing file into dataframe.

    Args:
        path: File path

    Returns: df
        DataFrame from table in file

    """
    sep = _check_suffix(path)
    df = pd.read_csv(_is_valid_file(path), sep=sep, index_col=0)
    return df


def parse_outliers(path: str, updown: str, iqrs: float) -> OutlierTable:
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
    outliers = OutlierTable(df, updown, iqrs, samples, None)
    outliers.make_frac_table()
    return outliers


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
                val = val.replace("_", "-")
                new_df.loc[(df[col] != val), binarized_col_name % (col, val)] = (
                    outgroup_val % val
                )
                new_df.loc[(df[col] == val), binarized_col_name % (col, val)] = val
                new_df.loc[df[col].isnull(), binarized_col_name % (col, val)] = np.nan
    return new_df


def list_to_file(lis: Iterable, filename: str):
    """Takes an iterable and a file path and writes a value per line from the iterable into the new
    file.

    Args:
        lis: Iterable to write to file
        filename: Filename to write to.

    Returns:
        None

    """

    with open(_check_output_prefix(filename), "w") as fh:
        for x in lis:
            fh.write("%s\n" % x)


def convert_to_residuals(
        valstarget: Series,
        valsnormer: Series,
        model
) -> Series:
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UndefinedMetricWarning)
    nonull = ((valstarget.isnull() == False) & (valsnormer.isnull() == False))
    if sum(nonull) < model.get_params()['cv']:
        residuals = np.empty(len(valstarget))
        residuals = pd.Series(residuals, index=valstarget.index)
    else:
        features = valsnormer[nonull].values.reshape(-1, 1)
        labels = valstarget[nonull].values
        model = model.fit(features, labels)
        prediction = model.predict(features)
        residuals = labels - prediction
        residuals = pd.Series(residuals, index=valstarget[nonull].index)
    return residuals


def normalize_df(
        target: DataFrame,
        normer: DataFrame,
        ind_sep: str = '-',
        alphas: Optional[Iterable[float]] = None,
        cv: float = 5,
        **RidgeCV_kws
) -> DataFrame:
    if alphas is None:
        alphas = [2 ** i for i in range(-10, 10, 1)]

    normer = normer.reindex(target.columns, axis=1)

    target = target.transpose()
    target['col0'] = 0
    target.set_index('col0', append=True, inplace=True)
    target = target.reorder_levels([target.index.names[-1], target.index.names[0]]).transpose()

    normer = normer.transpose()
    normer['col0'] = 1
    normer.set_index('col0', append=True, inplace=True)
    normer = normer.reorder_levels([normer.index.names[-1], normer.index.names[0]]).transpose()

    target['gene'] = [i.split(ind_sep)[0] for i in target.index]
    target = target.loc[target['gene'].isin(normer.index), :]
    data = target.merge(normer, how='left', left_on='gene', right_index=True)

    model = lm.RidgeCV(alphas=alphas, cv=cv, **RidgeCV_kws)
    normed = data.apply(
        (lambda row: convert_to_residuals(row[0], row[1], model)),
        axis=1)

    return normed