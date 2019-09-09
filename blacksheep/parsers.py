import warnings
import logging
import os.path
from typing import Iterable, Optional
import pandas as pd
from pandas import DataFrame, Series
import numpy as np
import sklearn.linear_model as lm
from sklearn.exceptions import UndefinedMetricWarning
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


def _convert_to_residuals(valstarget: Series, valsnormer: Series, model) -> Series:
    """ Take a 1d array of values to normalize, and a 1d array of values to use for
    normalization and a model. Trains the model using the target as labels and normer as
    features, then predicts for each normer values and subtracts the original target value,
    leaving residual values.

    Args:
        valstarget: 1d array of values to normalize
        valsnormer: 1d array of values to use to normalize
        model: a machine learning model, on which you can perform .fit and .predict

    Returns: residuals

    """
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=UndefinedMetricWarning)
    nonull = (valstarget.isnull() == False) & (valsnormer.isnull() == False)
    if sum(nonull) < model.get_params()["cv"]:
        residuals = np.empty(len(valstarget))
        return pd.Series(residuals, index=valstarget.index)
    features = valsnormer[nonull].values.reshape(-1, 1)
    labels = valstarget[nonull].values
    model = model.fit(features, labels)
    prediction = model.predict(features)
    residuals = labels - prediction
    return pd.Series(residuals, index=valstarget[nonull].index)


def normalize_df(
    target: DataFrame,
    normer: DataFrame,
    ind_sep: Optional[str] = "-",
    alphas: Optional[Iterable[float]] = None,
    cv: float = 5,
    **RidgeCV_kws
) -> DataFrame:
    """ Used to normalize a dataset by another dataset, using a linear model with regularization
    chosen through cross validation (aka sklearn's RidgeCV). This is useful for normalizing,
    for example, RNA values by CNA, or phosphopeptide values by protein abundance. If target and
    normer dataframe row IDs (index) match 1:1, pass None for ind_sep.

    Args:
        target: Dataframe of values to normalize. Row IDs (index) before the sep (or whole ID
        if no sep) must match normer IDs. Row IDs must be unique.
        normer: Dataframe of values to use for normalization. Row IDs must match all or
        pre-ind_sep portions of target row IDs. Row IDs must be unique.
        ind_sep: If multiple rows in target map to 1 row in normer, the delimiter used to split
        the unique ID that matches the normer IDs. Defaul "-"
        alphas: Parameters to try for regulariztion. If None, tries powers of 2 from -10 to 10.
        cv: Fold for cross validation. Also the minimum number of non-null values for each
        row. Default 5
        **RidgeCV_kws: kws to pass to sklearn's RidgeCV

    Returns: normed
        The target dataframe normalized by the normer dataframe. Only includes rows with
        sufficient non-null values from both dataframe.

    """

    if not alphas:
        alphas = [2 ** i for i in range(-10, 10, 1)]

    normer = normer[[col for col in target.columns if col in normer.columns]]
    target = target[normer.columns]
    if (len(normer.columns) < cv) or (len(target.columns) < cv):
        raise KeyError("target and normer dataframes do not have at least %s columns in common" %cv)

    target = target.transpose()
    target["col0"] = 0
    target.set_index("col0", append=True, inplace=True)
    target = target.reorder_levels(
        [target.index.names[-1], target.index.names[0]]
    ).transpose()

    normer = normer.transpose()
    normer["col0"] = 1
    normer.set_index("col0", append=True, inplace=True)
    normer = normer.reorder_levels(
        [normer.index.names[-1], normer.index.names[0]]
    ).transpose()

    target["gene"] = [i.split(ind_sep)[0] for i in target.index]
    target = target.loc[target["gene"].isin(normer.index), :]
    if len(target) == 0:
        raise KeyError("No rows in common between target and normer")
    logging.info(
        "Normalizing %s common rows and %s common samples between target and normer"
        %(len(target), len(normer.columns))
    )
    data = target.merge(normer, how="left", left_on="gene", right_index=True)

    model = lm.RidgeCV(alphas=alphas, cv=cv, **RidgeCV_kws)
    normed = data.apply(
        (lambda row: _convert_to_residuals(row[0], row[1], model)), axis=1
    )

    return normed
