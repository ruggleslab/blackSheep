from typing import Optional, List, Iterable
import sys
import os.path
import logging
import argparse
import matplotlib.pyplot as plt
from blacksheep.outliers import run_outliers
from blacksheep.outliers import make_outliers_table
from blacksheep.outliers import compare_groups_outliers
from blacksheep import parsers
from blacksheep.parsers import is_valid_file, check_output_prefix
from blacksheep.classes import qValues
from blacksheep.visualization import plot_heatmap
from blacksheep.constants import *


def set_up_logger(path):
    fmt = '%(asctime)s:%(levelname)s:%(message)s'
    logging.basicConfig(
        filename='test.log',
        format=fmt,
        level=logging.INFO,
        datefmt='%m/%d/%Y %H:%M:%S'
    )


# arg checkers
def check_positive(arg: str) -> float:
    try:
        arg = float(arg)
    except TypeError:
        raise argparse.ArgumentTypeError("%s is not a valid number" % arg)
    if arg < 0:
        raise argparse.ArgumentTypeError("%s is not a positive number" % arg)
    return arg


def bn0and1(arg: str) -> float:
    try:
        arg = float(arg)
    except TypeError:
        raise argparse.ArgumentTypeError("%s is an invalid positive value" % arg)
    if (arg < 0) or (arg > 1):
        raise argparse.ArgumentTypeError("%s is not a number between 0 and 1" % arg)
    return arg


# Argparser
def parse_args(args: List):
    parser = argparse.ArgumentParser(prog="BlackSheep", description="")
    parser.add_argument("--version", "-v", action="version", version="%(prog)s 0.0.1")

    subparsers = parser.add_subparsers(dest="which")
    subparsers.required = True

    outliers_table = subparsers.add_parser(
        "outliers_table",
        description="Takes a table of values and converts to a "
        "table of outlier counts.",
    )
    outliers_table.add_argument(
        "values",
        type=is_valid_file,
        help="File path to input values. First column must be sample names, "
        "header must be site labels. Only .tsv and .csv accepted.",
    )
    outliers_table.add_argument(
        "--iqrs",
        type=check_positive,
        default=1.5,
        help="Number of interquartile ranges (IQRs) above or below the "
        "median to consider a value an outlier. Default is 1.5 IQRs.",
    )
    outliers_table.add_argument(
        "--up_or_down",
        type=str,
        default="true",
        choices=["up", "down"],
        help="Whether to look for up or down outliers. Choices are up or "
        "down. Default up.",
    )
    outliers_table.add_argument(
        "--output_prefix",
        type=check_output_prefix,
        default="outliers_table",
        help="Output prefix for writing files. Default outliers. ",
    )
    outliers_table.add_argument(
        "--ind_sep",
        type=str,
        default="-",
        help="If site labels have a parent molecule (e.g. a gene name such "
        "as ATM) and a site identifier (e.g. S365) this is the "
        "delimiter between the two elements. Default is -",
    )
    outliers_table.add_argument(
        "--do_not_aggregate",
        default=False,
        action="store_true",
        help="Use flag if you do not want to sum outliers based on site " "prefixes.",
    )
    outliers_table.add_argument(
        "--write_frac_table",
        default=False,
        action="store_true",
        help="Use flag if you want to write a table with fraction of "
        "values per site, per sample that are outliers. Will not be "
        "written by default. Useful for visualization. ",
    )

    binarize = subparsers.add_parser(
        "binarize",
        description="Takes an annotation table where some columns "
        "may have more than 2 possible values (not "
        "including empty/null values) and outputs an "
        "annotation table with only two values per "
        "annotation. Propagates null values. ",
    )
    binarize.add_argument(
        "annotations",
        type=is_valid_file,
        help="Annotation table with "
        "samples as rows and "
        "annotation labels as "
        "columns. ",
    )
    binarize.add_argument(
        "--output_prefix",
        type=check_output_prefix,
        default="binarized_annotations",
        help="Prefix for new annotation file",
    )

    compare_groups = subparsers.add_parser(
        "compare_groups",
        description="Takes an annotation table and outlier count "
        "table (output of outliers_table) and outputs "
        "qvalues from a statistical test that looks "
        "for enrichment of outlier values in each "
        "group in the annotation table. For each value in each comparison, the qvalue table will "
        "have 1 column, if there are any genes in that comparison. ",
    )
    compare_groups.add_argument(
        "outliers_table",
        type=is_valid_file,
        help="Table of outlier counts (output of outliers_table). Must be "
        ".tsv or .csv file, with outlier and non-outlier counts as rows, and genes/sites as "
        "columns. ",
    )
    compare_groups.add_argument(
        "annotations",
        type=is_valid_file,
        help="Table of annotations. Must be .csv or .tsv. Samples as rows "
        "and comparisons as columns. Comparisons must have only  "
        "unique values (not including missing values). If there are "
        "more options than that, you can use binarize to prepare the "
        "table. ",
    )
    compare_groups.add_argument(
        "--frac_filter",
        type=bn0and1,
        default=0.3,
        help="The minimum fraction of samples per group that must have an outlier in a gene to"
        "consider that gene in the analysis. This is used to prevent a high number of outlier "
        "values in 1 sample from driving a low qvalue. Default 0.3",
    )
    compare_groups.add_argument(
        "--output_prefix",
        type=check_output_prefix,
        default="compare_groups",
        help="File prefix for qvalues table and other optional outputs. Default outliers",
    )
    compare_groups.add_argument(
        "--write_comparison_summaries", default=False, action="store_true",
        help="Use flag to write a separate file for each column in the annotations table, "
             "with outlier counts in each group, p-values and q-values in each group. "
    )
    compare_groups.add_argument(
        "--iqrs", type=check_positive, default=None,
        help="Number of IQRs used to define outliers in the input count table. Optional."
    )
    compare_groups.add_argument(
        "--up_or_down", type=str, choices=["up", "down"],
        help="Whether input outlier table represents up or down outliers. Needed for "
             "output file labels. Default up"
    )
    compare_groups.add_argument(
        "--write_gene_list", default=False, action="store_true",
        help="Use flag to write a list of significantly enriched genes for each value in each "
             "comparison. If used, need an fdr threshold as well. "
    )
    compare_groups.add_argument(
        "--make_heatmaps", default=False, action="store_true",
        help="Use flag to draw a heatmap of signficantly enriched genes for each value in each "
             "comparison. If used, need an fdr threshold as well. "
    )
    compare_groups.add_argument(
        "--fdr", type=bn0and1, default=0.05,
        help="FDR cut off to use for signficantly enriched gene lists and heatmaps. Default 0.05"
    )
    compare_groups.add_argument(
        "--red_or_blue", type=str, choices=["red", "blue"], default="red",
        help="If --make_heatmaps is called, color of values to draw on heatmap. Default red. "
    )
    compare_groups.add_argument(
        "--annotation_colors", type=is_valid_file, default=None,
        help="File with color map to use for annotation header if --make_heatmaps is used. Must "
             "have a 'value    color' format for each value in annotations. Any value not "
             "represented will be assigned a new color. "
    )

    visualize = subparsers.add_parser("visualize",
                                      description="Used to make custom heatmaps from significant "
                                                "genes. ")
    visualize.add_argument(
        "comparison_qvalues", type=is_valid_file,
        help="Table of qvalues, output from compare_groups. Must be .csv or .tsv. Has genes/sites "
             "as rows and comparison values as columns. "
    )
    visualize.add_argument(
        "annotations", type=is_valid_file, help="Table of annotations used to generate qvalues. "
    )
    visualize.add_argument(
        "visualization_table", type=is_valid_file,
        help="Values to visualize in heatmap. Samples as rows and "
             "genes/sites as columns. Using outlier fraction table is recommended, but original "
             "values can also be used if no aggregation was used. "
    )
    visualize.add_argument(
        "comparison_of_interest", type=str,
        help="Name of column in qvalues table from which to visualize significant genes. "
    )
    visualize.add_argument(
        "--annotations_to_show", type=str, default=None, nargs="+",
        help="Names of columns from the annotation table to show in the header of the heatmap. "
             "Default is all columns. "
    )
    visualize.add_argument("--fdr", type=bn0and1, default=0.05,
                           help="FDR threshold to use to select genes to visualize. Default 0.05")
    visualize.add_argument(
        "--red_or_blue", type=str, choices=["red", "blue"], default="red",
        help="Color of values to draw on heatmap. Default red. "
    )
    visualize.add_argument(
        "--output_prefix", type=check_output_prefix, default="visualization",
        help="Output prefix for writing files. Default outliers. "
    )
    visualize.add_argument(
        "--annotation_colors", type=is_valid_file, default=None,
        help="File with color map to use for annotation header. Must "
             "have a line with 'value    color' format for each value in annotations. Any value "
             "not represented will be assigned a new color. "
    )
    visualize.add_argument(
        "--write_gene_list", default=False, action="store_true",
        help="Use flag to write a list of significantly enriched genes for each value in each "
             "comparison."
    )

    outliers = subparsers.add_parser(
        "outliers", description="Runs whole outliers pipeline. Has options to output every "
                                "possible output. "
    )
    outliers.add_argument(
        "values",
        type=is_valid_file,
        help="File path to input values. Rows are sample names, "
        "columns are gene/site labels. Only .tsv and .csv accepted.",
    )
    outliers.add_argument(
        "annotations",
        type=is_valid_file,
        help="File path to annotation values. Rows are sample names, "
             "header is different annotations. e.g. mutation status.",
    )
    outliers.add_argument(
        "--iqrs",
        type=check_positive,
        default=1.5,
        help="Number of inter-quartile ranges (IQRs) above or below the "
             "median to consider a value an outlier. Default is 1.5.",
    )
    outliers.add_argument(
        "--up_or_down",
        type=str,
        default="true",
        choices=["up", "down"],
        help="Whether to look for up or down outliers. Choices are up or "
        "down. Default up.",
    )
    outliers.add_argument(
        "--do_not_aggregate",
        default=False,
        action="store_true",
        help="Use flag if you do not want to sum outliers based on site prefixes."
    )
    outliers.add_argument(
        "--output_prefix",
        type=check_output_prefix,
        default="outliers_pipeline",
        help="Output prefix for writing files. Default outliers. ",
    )
    outliers.add_argument(
        "--write_outlier_table", default=False, action="store_true",
        help="Use flag to write a table of outlier counts."
    )
    outliers.add_argument(
        "--write_frac_table",
        default=False,
        action="store_true",
        help="Use flag if you want to write a table with fraction of "
        "values per site per sample that are outliers. Useful for custom visualization."
    )
    outliers.add_argument(
        "--ind_sep",
        type=str,
        default="-",
        help="If site labels have a parent molecule (e.g. a gene name such "
        "as ATM) and a site identifier (e.g. S365) this is the "
        "delimiter between the two elements. Default is -"
    )
    outliers.add_argument(
        "--frac_filter", type=bn0and1, default=0.3,
        help="The minimum fraction of samples per group that must have an outlier in a gene to"
             "consider that gene in the analysis. This is used to prevent a high number of outlier "
             "values in 1 sample from driving a low qvalue. Default 0.3"
    )
    outliers.add_argument(
        "--write_comparison_summaries", default=False, action="store_true",
        help="Use flag to write a separate file for each column in the annotations table, "
             "with outlier counts in each group, p-values and q-values in each group. "
    )
    outliers.add_argument("--fdr", type=bn0and1, default=0.05,
                           help="FDR threshold to use to select genes to visualize. Default 0.05")
    outliers.add_argument(
        "--write_gene_list", default=False, action="store_true",
        help="Use flag to write a list of significantly enriched genes for each value in each "
             "comparison."
    )
    outliers.add_argument(
        "--make_heatmaps", default=False, action="store_true",
        help="Use flag to draw a heatmap of significantly enriched genes for each value in each "
             "comparison. If used, need an fdr threshold as well. "
    )
    outliers.add_argument(
        "--red_or_blue", type=str, choices=["red", "blue"], default="red",
        help="Color of values to draw on heatmap. Default red. "
    )
    outliers.add_argument(
        "--annotation_colors", type=is_valid_file, default=None,
        help="File with color map to use for annotation header. Must "
             "have a line with 'value    color' format for each value in annotations. Any value "
             "not represented will be assigned a new color. "
    )
    return parser.parse_args(args)


# Module runner
def main(args: Optional[List[str]] = None):
    if args is None:
        args = sys.argv[1:]
    args = parse_args(args)

    set_up_logger(args.output_prefix)
    logging.info('Running BlackSheep in %s mode' % args.which)
    for arg in vars(args):
        logging.info("Parameter\n%s: %s" % (arg, getattr(args, arg)))

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
        annotations = parsers.parse_values(args.annotations)
        annotations = parsers.binarize_annotations(annotations)
        annotations.to_csv("%s.tsv" % args.output_prefix, sep="\t")

    elif args.which == "compare_groups":
        outliers = parsers.parse_outliers(
            args.outliers_table, args.up_or_down, args.iqrs
        )

        annotations = parsers.parse_values(args.annotations)
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
        annotations = parsers.parse_values(args.annotations)
        frac_table = parsers.parse_values(args.visualization_table)
        col_of_interest = args.comparison_of_interest
        annot_cols = args.annotations_to_show[0].split()

        try:
            annotations = annotations[annot_cols]
        except KeyError:
            raise ValueError("Some of %s not in annotations" % annot_cols)

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
            qValues(qvals[[col_of_interest]], [], None).write_gene_lists(
                args.fdr, args.output_prefix
            )

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

    with open(parameters_file_name % args.output_prefix, "w") as fh:
        for arg in vars(args):
            fh.write("%s: %s\n" % (arg, getattr(args, arg)))


if __name__ == "__main__":
    main()
