import argparse
from blacksheep.cli import main

def test_cli_outliers_table():
    pass


def test_cli_compare_groups():
    pass


def test_cli_outliers():
    parser = main()
    parser.parse_args(['-g', 'xyz', 'foo', '--count', '42'])
    pass
