import pickle
import blacksheep as bsh


def test_outliers_table():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)

    test_outliers = bsh.make_outliers_table(df)
    test_df = test_outliers.df.sort_index().sort_index(axis=1)
    outliers = outliers.sort_index().sort_index(axis=1)
    assert outliers.equals(test_df)


def test_frac_table():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)

    test_outliers = bsh.make_outliers_table(df)
    test_frac = test_outliers.frac_table

    test_frac = test_frac.sort_index().sort_index(axis=1)
    fractable = fractable.sort_index().sort_index(axis=1)
    assert test_frac.equals(fractable)


def test_compare_groups_qvals():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)
    outliers = bsh.classes.OutlierTable(outliers, "up", 1.5, df.columns, fractable)
    test_qvals = bsh.compare_groups_outliers(outliers, annotations)
    test_qvals.df = test_qvals.df.sort_index().sort_index(axis=1)
    qvalues = qvalues.sort_index().sort_index(axis=1)
    assert qvalues.equals(test_qvals.df)


def test_compare_groups_comps():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)

    outliers = bsh.classes.OutlierTable(outliers, "up", 1.5, df.columns, fractable)
    test_qvals = bsh.compare_groups_outliers(outliers, annotations)

    assert sum(annotations.columns.sort_values() != test_qvals.comps.sort_values()) == 0


def test_run_outliers_outTable():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)

    test_outliers, test_qvals = bsh.deva(df, annotations)
    outliers = outliers.sort_index().sort_index(axis=1)
    test_outliers.df = test_outliers.df.sort_index().sort_index(axis=1)

    assert outliers.equals(test_outliers.df)


def test_run_outliers_fracTable():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)

    test_outliers, test_qvals = bsh.deva(df, annotations)
    fractable = fractable.sort_index().sort_index(axis=1)
    test_outliers.frac_table = test_outliers.frac_table.sort_index().sort_index(axis=1)
    assert fractable.equals(test_outliers.frac_table)


def test_run_outliers_qvals():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)

    _, test_qvals = bsh.deva(df, annotations)
    test_qvals.df = test_qvals.df.sort_index().sort_index(axis=1)
    qvalues = qvalues.sort_index().sort_index(axis=1)

    assert qvalues.equals(test_qvals.df)


def test_run_outliers_comps():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)

    test_outliers, test_qvals = bsh.deva(df, annotations)
    assert sum(annotations.columns.sort_values() != test_qvals.comps.sort_values()) == 0


def test_normalizer():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, fractable, qvalues = pickle.load(fh)
    bsh.normalize(df)

