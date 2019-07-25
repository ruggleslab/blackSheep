import os
import pandas as pd
from blacksheep import visualization


def testHeatmap():
    annotations = pd.read_csv('tests/sample_annotations.csv', index_col=0)
    qvals = pd.read_csv('tests/output/pipeline_test.up.qvalues.tsv', sep='\t', index_col=0)
    col_of_interest = 'fisherFDR_Histologic_type_Serous'
    fdr = 0.5
    vis_table = pd.read_csv('tests/output/pipeline_test.up.fraction_table.tsv', sep='\t',
                            index_col=0)
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


