import sys
import os

sys.path.append("..")
import pandas as pd
from blacksheep import *


def test_make_outliers_table():

    input_df = pd.read_csv("sample_endo.csv", index_col=0)
    outliers = make_outliers_table(input_df)
    test_df = outliers.df

    standard_df = pd.read_csv("sample_endo_outliers_table.csv", index_col=0)

    assert test_df.equals(standard_df)
