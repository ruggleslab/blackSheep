import pandas as pd
import numpy as np
import blacksheep as bsh


def test_make_outliers_table():

    input_df = pd.read_csv("sample_endo.csv", index_col=0)
    outliers = bsh.make_outliers_table(input_df)
    test_df = outliers.df

    standard_df = pd.read_csv("sample_endo_outliers_table.csv", index_col=0)

    assert test_df.equals(standard_df)


def test_make_outliers_table():
    values = np.array[[np.nan, -2, -1, -0.5, 0.1, 0.25, 0.35, 1, 4.5]]