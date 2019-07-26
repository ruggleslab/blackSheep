from typing import Optional, List
import sys
import logging
import argparse
import os.path
import matplotlib.pyplot as plt
from .outliers import run_outliers
from .outliers import make_outliers_table
from .outliers import compare_groups_outliers
from . import parsers
from .visualization import plot_heatmap
from .visualization import write_genes


# File checkers
def check_positive(arg):
    try:
        arg = float(arg)
    except TypeError:
        raise argparse.ArgumentTypeError("%s is not a valid number" % arg)
    if arg < 0:
        raise argparse.ArgumentTypeError("%s is not a positive number" % arg)
    return arg


def is_valid_file(arg):
    if not os.path.exists(arg):
        raise argparse.ArgumentTypeError("%s does not exist" % arg)
    return arg


def check_output_prefix(arg):
    if "/" in arg:
        prefix = arg.rsplit("/", 1)[0]
        is_valid_file(prefix)
    return arg


def bn0and1(arg):
    try:
        arg = float(arg)
    except TypeError:
        raise argparse.ArgumentTypeError("%s is an invalid positive value" % arg)
    if (arg < 0) or (arg > 1):
        raise argparse.ArgumentTypeError("%s is not a number between 0 and 1" % arg)
    return arg


def cols_file_or_list(arg):
    arg = arg.split()

    if len(arg) > 1:
        return arg
    if len(arg) != 1:
        raise argparse.ArgumentTypeError("Cannot read %s as list or file" % arg)
    try:
        with open(arg[0], "r") as fh:
            return [x for x in fh.readlines()]
    except TypeError:
        raise argparse.ArgumentTypeError("Cannot read %s as list or file" % arg)

# Argparser
def parse_args(args):
    # TODO finish helps
    parser = argparse.ArgumentParser(prog="BlackSheep")
    parser.add_argument("--version", action="version", version="%(prog)s 0.0.1")

    subparsers = parser.add_subparsers(dest="which")
    subparsers.required = True

    mk_out_table = subparsers.add_parser("outliers_table")
    mk_out_table.add_argument(
        "values",
        type=is_valid_file,
        help="File path to input values. First column must be sample names, "
        "header must be site labels. Only .tsv and .csv accepted.",
    )
    mk_out_table.add_argument(
        "--iqrs",
        type=check_positive,
        default=1.5,
        help="Number of interquartile ranges (IQRs) above or below the "
        "median to consider a value an outlier. Default is 1.5 IQRs.",
    )
    mk_out_table.add_argument(
        "--up_or_down",
        type=str,
        default="true",
        choices=["up", "down"],
        help="Whether to look for up or down outliers. Choices are up or "
        "down. Default up.",
    )
    mk_out_table.add_argument(
        "--output_prefix",
        type=check_output_prefix,
        default="outliers",
        help="Output prefix for writing files. Default outliers. ",
    )
    mk_out_table.add_argument(
        "--ind_sep",
        type=str,
        default="-",
        help="If site labels have a parent molecule (e.g. a gene name such "
        "as ATM) and a site identifier (e.g. S365) this is the "
        "delimiter between the two elements. Default is -",
    )
    mk_out_table.add_argument(
        "--do_not_aggregate",
        default=False,
        action="store_true",
        help="Use flag if you do not want to sum outliers based on site " "prefixes.",
    )
    mk_out_table.add_argument(
        "--write_frac_table",
        default=False,
        action="store_true",
        help="Use flag if you want to write a table with fraction of "
        "values per site, per sample that are outliers. Will not be "
        "written by default. Useful for visualization. ",
    )

    binarize = subparsers.add_parser("binarize")
    binarize.add_argument("annotation_table", type=is_valid_file)
    binarize.add_argument("--new_file_prefix", type=check_output_prefix,
                          default="binarized_annotations")

    comparisons = subparsers.add_parser("compare_groups")
    comparisons.add_argument("outliers_table", type=is_valid_file)
    comparisons.add_argument("annotations_table", type=is_valid_file)
    comparisons.add_argument("--frac_filter", type=bn0and1, default=0.3)
    comparisons.add_argument(
        "--output_prefix", type=check_output_prefix, default="outliers"
    )
    comparisons.add_argument(
        "--write_comparison_summaries", default=False, action="store_true"
    )
    comparisons.add_argument(
        "--iqrs", type=check_positive, default=1.5, help="Number of IQRs in "
    )
    comparisons.add_argument("--up_or_down", type=str, choices=["up", "down"])
    comparisons.add_argument("--fdr", type=bn0and1, default=0.05)
    comparisons.add_argument("--write_gene_list", default=False, action="store_true")
    comparisons.add_argument("--make_heatmaps", default=False, action="store_true")
    comparisons.add_argument(
        "--red_or_blue", type=str, choices=["red", "blue"], default="red"
    )
    comparisons.add_argument("--annotation_colors", type=is_valid_file, default=None)

    figures = subparsers.add_parser("visualize")
    figures.add_argument("comparison_qvalues", type=is_valid_file)
    figures.add_argument("annotation_header", type=is_valid_file)
    figures.add_argument("fraction_table", type=is_valid_file)
    figures.add_argument("comparison_of_interest", type=str)
    figures.add_argument(
        "--annotations_to_show", type=cols_file_or_list, default=None, nargs="+"
    )
    figures.add_argument("--fdr", type=bn0and1, default=0.05)
    figures.add_argument(
        "--red_or_blue", type=str, choices=["red", "blue"], default="red"
    )
    figures.add_argument(
        "--output_prefix", type=check_output_prefix, default="outliers"
    )
    figures.add_argument("--annotation_colors", type=is_valid_file, default=None)
    figures.add_argument("--write_gene_list", default=False, action="store_true")

    pipeline = subparsers.add_parser("outliers")
    pipeline.add_argument(
        "values",
        type=is_valid_file,
        help="File path to input values. First column must be sample names, "
        "header must be site labels. Only .tsv and .csv accepted.",
    )
    pipeline.add_argument(
        "annotations",
        type=is_valid_file,
        help="File path to annotation values. First column must be sample "
        "names, header is different annotations. e.g. mutation status.",
    )
    pipeline.add_argument(
        "--iqrs",
        type=check_positive,
        default=1.5,
        help="Number of inter-quartile ranges (IQRs) above or below the "
        "median to consider a value an outlier. Default is 1.5 IQRs.",
    )
    pipeline.add_argument(
        "--up_or_down",
        type=str,
        default="true",
        choices=["up", "down"],
        help="Whether to look for up or down outliers. Choices are up or "
        "down. Default up.",
    )
    pipeline.add_argument(
        "--do_not_aggregate",
        default=False,
        action="store_true",
        help="Use flag if you do not want to sum outliers based on site " "prefixes.",
    )
    pipeline.add_argument(
        "--output_prefix",
        type=check_output_prefix,
        default="outliers",
        help="Output prefix for writing files. Default outliers. ",
    )
    pipeline.add_argument(
        "--write_frac_table",
        default=False,
        action="store_true",
        help="Use flag if you want to write a table with fraction of "
        "values per site, per sample that are outliers. Will not be "
        "written by default. Useful for visualization. ",
    )
    pipeline.add_argument(
        "--write_outlier_table", default=False, action="store_true", help=""
    )
    pipeline.add_argument(
        "--ind_sep",
        type=str,
        default="-",
        help="If site labels have a parent molecule (e.g. a gene name such "
        "as ATM) and a site identifier (e.g. S365) this is the "
        "delimiter between the two elements. Default is -",
    )
    pipeline.add_argument("--frac_filter", type=bn0and1, default=0.3)
    pipeline.add_argument(
        "--write_comparison_summaries", default=False, action="store_true"
    )
    pipeline.add_argument("--fdr", type=bn0and1, default=0.05)
    pipeline.add_argument("--write_gene_list", default=False, action="store_true")
    pipeline.add_argument("--make_heatmaps", default=False, action="store_true")
    pipeline.add_argument(
        "--red_or_blue", type=str, choices=["red", "blue"], default="red"
    )
    pipeline.add_argument("--annotation_colors", type=is_valid_file, default=None)
    return parser.parse_args(args)


# Module runner
def main(args: Optional[List[str]] = None):
    if args is None:
        args = sys.argv[1:]
    args = parse_args(args)

    if args.which == "outliers_table":
        df = parsers.parse_values(args.values)
        make_outliers_table(
            df,
            args.iqrs,
            args.up_or_down,
            ~args.do_not_aggregate,
            save_outlier_table=True,
            save_frac_table=args.write_frac_table,
            output_prefix=args.output_prefix,
            ind_sep=args.ind_sep,
        )

    elif args.which == "binarize":
        annotations = parsers.parse_values(args.annotation_table)
        annotations = parsers.binarize_annotations(annotations)
        annotations.to_csv('%s.tsv' % args.output, sep='\t')

    elif args.which == "compare_groups":
        outliers = parsers.parse_outliers(
            args.outliers_table, args.up_or_down, args.iqrs
        )

        annotations = parsers.parse_values(args.annotations_table)
        qVals = compare_groups_outliers(
            outliers,
            annotations,
            args.frac_filter,
            save_qvalues=True,
            output_prefix=args.output_prefix,
            output_comparison_summaries=args.write_comparison_summaries,
        )
        if args.write_gene_list:
            qVals.write_gene_lists(args.fdr, args.output_prefix)

        if args.make_heatmaps:
            for col_of_interest in qVals.df.columns:
                plot_heatmap(
                    annotations,
                    qVals.df,
                    col_of_interest,
                    outliers.frac_table,
                    fdr=args.fdr,
                    red_or_blue=args.red_or_blue,
                    output_prefix=args.output_prefix,
                    colors=args.annotation_colors,
                    savefig=True,
                )
                plt.close()

    elif args.which == "visualize":
        qvals = parsers.parse_values(args.comparison_qvalues)
        annotations = parsers.parse_values(args.annotation_header)
        frac_table = parsers.parse_values(args.fraction_table)
        col_of_interest = args.comparison_of_interest
        annot_cols = args.annotations_to_show[0]

        try:
            annotations = annotations[annot_cols]
        except KeyError:
            print(annotations.columns, annotations.index, annotations.head())
            raise ValueError("Some cols not in annotations: %s" % annot_cols)

        annot_label = col_of_interest.split("_", 1)[1].rsplit("_", 1)[0]
        if annot_label not in annot_cols:
            raise ValueError("%s must be in annotations columns to show" % annot_label)

        plot_heatmap(
            annotations,
            qvals,
            col_of_interest,
            frac_table,
            fdr=args.fdr,
            red_or_blue=args.red_or_blue,
            output_prefix=args.output_prefix,
            colors=args.annotation_colors,
            savefig=True,
        )
        plt.close()
        if args.write_gene_list:
            write_genes(qvals[[col_of_interest]], args.fdr, args.output_prefix)

    elif args.which == "outliers":
        df = parsers.parse_values(args.values)
        annotations = parsers.parse_values(args.annotations)
        outLiers, qVals = run_outliers(
            df,
            annotations,
            args.frac_filter,
            args.iqrs,
            args.up_or_down,
            ~args.do_not_aggregate,
            args.write_outlier_table,
            args.write_frac_table,
            save_qvalues=True,
            output_prefix=args.output_prefix,
            ind_sep=args.ind_sep,
            output_comparison_summaries=args.write_comparison_summaries,
        )
        if args.write_gene_list:
            qVals.write_gene_lists(args.fdr, args.output_prefix)

        if args.make_heatmaps:
            for col_of_interest in qVals.df.columns:
                plot_heatmap(
                    annotations,
                    qVals.df,
                    col_of_interest,
                    outLiers.frac_table,
                    fdr=args.fdr,
                    red_or_blue=args.red_or_blue,
                    output_prefix=args.output_prefix,
                    colors=args.annotation_colors,
                    savefig=True,
                )
                plt.close()

    with open("%s.parameters.txt" % args.output_prefix, "w") as fh:
        for arg in vars(args):
            fh.write("%s: %s\n" % (arg, getattr(args, arg)))


if __name__ == "__main__":
    main()
