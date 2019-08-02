import pickle
import pandas as pd
from blacksheep import visualization


def test_heatmap_small():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, vis_table, qvals = pickle.load(fh)
    col_of_interest = "fisherFDR_comp0_1"
    fdr = 0.5

    red_or_blue = "red"
    output_prefix = "tests/output/vis_test_small"
    colors = None

    visualization.plot_heatmap(
        annotations,
        qvals,
        col_of_interest,
        vis_table,
        fdr,
        red_or_blue,
        output_prefix,
        colors,
        savefig=True,
    )


def test_comps_heatmap_large():
    annotations = pd.read_csv("tests/sample_annotations.csv", index_col=0)
    qvals = pd.read_csv("tests/sample_endo.up.qvalues.tsv", sep="\t", index_col=0)
    vis_table = pd.read_csv(
        "tests/sample_endo.up.fraction_table.tsv", sep="\t", index_col=0
    )
    col_of_interest = "fisherFDR_POLE_subtype_Yes"
    red_or_blue = "blue"
    output_prefix = "tests/output/vis_test_large"
    colors = None
    for fdr in [0.6, 0.05]:
        visualization.plot_heatmap(
            annotations,
            qvals,
            col_of_interest,
            vis_table,
            fdr,
            red_or_blue,
            output_prefix,
            colors,
            savefig=True,
        )


def test_comps_heatmap_medium():
    annotations = pd.read_csv("tests/sample_annotations.csv", index_col=0)
    qvals = pd.read_csv("tests/sample_endo.up.qvalues.tsv", sep="\t", index_col=0)
    vis_table = pd.read_csv(
        "tests/sample_endo.up.fraction_table.tsv", sep="\t", index_col=0
    )
    col_of_interest = "fisherFDR_MSI_status_MSI-H"
    red_or_blue = "blue"
    output_prefix = "tests/output/vis_test_medium"
    colors = None
    for fdr in [0.6, 0.05]:
        visualization.plot_heatmap(
            annotations,
            qvals,
            col_of_interest,
            vis_table,
            fdr,
            red_or_blue,
            output_prefix,
            colors,
            savefig=True,
        )
