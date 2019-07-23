import pandas as pd
from blacksheep.cli import main

def test_cli_outliers_table():
    args = ['outliers_table', 'tests/sample_endo.csv', '--iqrs', '1.5',
            '--up_or_down', 'up',
            '--output_prefix', 'tests/output/outliers_table_test',
            '--ind_sep', '-', '--write_frac_table']

    main(args)


def test_cli_compare_groups():
    args = ['compare_groups', 'tests/sample_endo_outliers_table.csv',
            'tests/sample_annotations.csv', '--iqrs', '1.5',
            '--up_or_down', 'up',
            '--output_prefix', 'tests/output/compare_groups_test',
            '--frac_filter', '0.1', '--write_comparison_summaries']

    main(args)


def test_cli_pipeline():
    args = ['outliers', 'tests/sample_endo.csv', 'tests/sample_annotations.csv', '--iqrs', '1.5',
            '--up_or_down', 'up',
            '--output_prefix', 'tests/output/pipeline_test', '--write_frac_table',
            '--write_outlier_table',
            '--ind_sep', '-', '--frac_filter', '0.10', '--write_comparison_summaries']

    main(args)

    # stnd_frac_table = pd.read_csv('')
    # stnd_outlier_table = pd.read_csv('sample_endo_outliers_table.csv', index_col=0)
    # stnd_qvals = pd.read_csv()
    pass
