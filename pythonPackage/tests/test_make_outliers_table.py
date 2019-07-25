import pandas as pd
import blacksheep as bsh

def test_make_outliers_table():

    input_df = pd.read_csv("sample_endo.csv", index_col=0)
    outliers = bsh.make_outliers_table(input_df)
    test_df = outliers.df
    # test_frac_table = outliers.frac_table
    standard_df = pd.read_csv("sample_endo_outliers_table.csv", index_col=0)

    assert test_df.equals(standard_df)


def test_make_outliers_table():

    pass