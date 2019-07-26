import pickle
import pandas as pd
from blacksheep import visualization


def test_heatmap_small():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, vis_table, qvals = pickle.load(fh)
    col_of_interest = "fisherFDR_comp0_1"
    fdr = 0.5

    red_or_blue = "red"
    output_prefix = "tests/output/vis_test"
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
    )


def test_comps_heatmap_larger():
    annotations = pd.read_csv("tests/sample_annotations.csv", index_col=0)
    qvals = pd.read_csv("tests/sample_endo.up.qvalues.tsv", sep="\t", index_col=0)
    vis_table = pd.read_csv(
        "tests/sample_endo.up.fraction_table.tsv", sep="\t", index_col=0
    )
    col_of_interest = "fisherFDR_Histologic_type_Serous"
    red_or_blue = "red"
    output_prefix = "tests/output/vis_test"
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
        )
