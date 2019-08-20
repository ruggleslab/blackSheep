from blacksheep.cli import _main


def test_cli_outliers_table():
    args = [
        "outliers_table",
        "tests/pidgin_values.csv",
        "--iqrs",
        "1.5",
        "--up_or_down",
        "up",
        "--output_prefix",
        "tests/output/outliers_table_test",
        "--ind_sep",
        "-",
        "--write_frac_table",
    ]

    _main(args)


def test_cli_normer():
    args = [
        "normalize",
        "tests/pidgin_values.csv",
        "tests/pidgin_fracTable.csv",
        "--output_prefix",
        "tests/output/normer",
        "--ind_sep",
        "-",
    ]

    _main(args)


def test_cli_compare_groups_with_heatmaps():
    args = [
        "compare_groups",
        "tests/pidgin_outliers.csv",
        "tests/pidgin_annotations.csv",
        "--iqrs",
        "1.5",
        "--up_or_down",
        "up",
        "--output_prefix",
        "tests/output/compare_groups_test",
        "--frac_filter",
        "0.1",
        "--write_comparison_summaries",
        "--make_heatmaps",
        "--write_gene_list",
        "--fdr",
        "0.6",
        "--red_or_blue",
        "red",
    ]

    _main(args)


def test_vis():
    args = [
        "visualize",
        "tests/pidgin_qvalues.csv",
        "tests/pidgin_annotations.csv",
        "tests/pidgin_fracTable.csv",
        "fisherFDR_comp0_1",
        "--annotations_to_show",
        "comp0 comp1",
        "--output_prefix",
        "tests/output/vis_cli_test",
        "--write_gene_list",
        "--fdr",
        "0.5",
        "--red_or_blue",
        "red",
    ]

    _main(args)


def test_cli_pipeline():
    args = [
        "deva",
        "tests/pidgin_values.csv",
        "tests/pidgin_annotations.csv",
        "--iqrs",
        "1.5",
        "--up_or_down",
        "up",
        "--output_prefix",
        "tests/output/pipeline_test",
        "--write_frac_table",
        "--write_outlier_table",
        "--ind_sep",
        "-",
        "--frac_filter",
        "0.10",
        "--write_comparison_summaries",
    ]

    _main(args)


def test_cli_pipeline_with_figures():
    args = [
        "deva",
        "tests/pidgin_values.csv",
        "tests/pidgin_annotations.csv",
        "--iqrs",
        "1.5",
        "--up_or_down",
        "up",
        "--output_prefix",
        "tests/output/pipeline_figs_test",
        "--ind_sep",
        "-",
        "--frac_filter",
        "0.10",
        "--write_comparison_summaries",
        "--fdr",
        "0.5",
        "--write_gene_list",
        "--make_heatmaps",
        "--red_or_blue",
        "blue",
    ]

    _main(args)


