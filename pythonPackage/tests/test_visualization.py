import pickle
from blacksheep import visualization


def test_heatmap():
    with open("tests/pidgin_example.pickle", "rb") as fh:
        df, annotations, outliers, vis_table, qvals = pickle.load(fh)
    col_of_interest = 'fisherFDR_comp0_1'
    fdr = 0.5

    red_or_blue = 'red'
    output_prefix = 'tests/output/vis_test'
    colors = None

    visualization.plot_heatmap(
        annotations,
        qvals,
        col_of_interest,
        vis_table,
        fdr,
        red_or_blue,
        output_prefix,
        colors
    )


