import pandas as pd
import numpy as np
from .classes import OutlierTable
from .outlierTable import makeFracTable


def parseValues(path):
    if path[-3:] == "tsv":
        sep = "\t"
    elif path[-3:] == "csv":
        sep = ","
    else:
        raise ValueError("File must be .csv or .tsv")
    df = pd.read_csv(path, sep=sep, index_col=0)
    return df


def parseOutliers(path, updown, iqrs):
    if path[-3:] == "tsv":
        sep = "\t"
    elif path[-3:] == "csv":
        sep = ","
    else:
        raise ValueError("File must be .csv or .tsv")
    df = pd.read_csv(path, sep=sep, index_col=0)
    samples = sorted(list(set([ind.rsplit("_", 1)[0] for ind in df.index])))
    frac_table = makeFracTable(df, samples)
    return OutlierTable(df, updown, iqrs, samples, frac_table)


def parseAnnotations(path):
    if path[-3:] == "tsv":
        sep = "\t"
    elif path[-3:] == "csv":
        sep = ","
    else:
        raise ValueError("File must be .csv or .tsv")
    df = pd.read_csv(path, sep=sep, index_col=0)
    return df


def binarizeAnnotations(df):
    new_df = pd.DataFrame(index=df.index)
    for col in df.columns:
        if len(df[col].dropna().value_counts().keys()) == 2:
            new_df[col] = df[col]
        elif len(df[col].dropna().value_counts().keys()) > 2:
            for val in df[col].dropna().value_counts().keys():
                new_df.loc[(df[col] != val), "%s_%s" % (col, val)] = "not_%s" % val
                new_df.loc[(df[col] == val), "%s_%s" % (col, val)] = val
                new_df.loc[df[col].isnull(), "%s_%s" % (col, val)] = np.nan
    return new_df
