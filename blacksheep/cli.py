from typing import Optional, List, Iterable
import sys
import logging
import argparse

import matplotlib.pyplot as plt
from blacksheep.deva import deva
from blacksheep.deva import make_outliers_table
from blacksheep.deva import compare_groups_outliers
from blacksheep import parsers
from blacksheep.parsers import _is_valid_file, _check_output_prefix, subset_by_genes
from blacksheep.classes import qValues
from blacksheep.simulate import run_simulations
from blacksheep.visualization import plot_heatmap
from blacksheep._constants import *

fmt = "%(asctime)s:%(levelname)s:%(message)s"
logging.basicConfig(format=fmt, level=logging.INFO, datefmt="%m/%d/%Y %H:%M:%S")
logging.captureWarnings(True)

def _set_up_logger(path):
    logger = logging.getLogger("cli")
    fh = logging.FileHandler("%s.log" % path)
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    fmter = logging.Formatter(fmt)
    fh.setFormatter(fmter)
    ch.setFormatter(fmter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# arg checkers
def _check_positive(arg: str) -> float:
    try:
        arg = float(arg)
    except TypeError:
        raise argparse.ArgumentTypeError("%s is not a valid number" % arg)
    if arg < 0:
        raise argparse.ArgumentTypeError("%s is not a positive number" % arg)
    return arg


def _bn0and1(arg: str) -> float:
    try:
        arg = float(arg)
    except TypeError:
        raise argparse.ArgumentTypeError("%s is an invalid positive value" % arg)
    if (arg < 0) or (arg > 1):
        raise argparse.ArgumentTypeError("%s is not a number between 0 and 1" % arg)
    return arg


# Argparser
def _make_parser():
    parser = argparse.ArgumentParser(prog="blacksheep", description="")
    parser.add_argument("--version", "-v", action="version", version="%(prog)s 0.0.1")

    subparsers = parser.add_subparsers(dest="which")
    subparsers.required = True

    normalize = subparsers.add_parser(
        "normalize",
        description="Takes an unnormalized values table and uses median of ratios normalization "
                    "to normalize. Saves a log2 normalized table appropriate for BlackSheep "
                    "analysis."
    )
    normalize.add_argument(
        "unnormed_values",
        type=_is_valid_file,
        help="Table of values to be normalized. Sites/genes as rows, samples as columns. "
    )
    normalize.add_argument(
        "--output_prefix",
        type=str,
        default="values",
        help="Prefix for output file. Suffix will be '.normalized.tsv'"
    )

    outliers_table = subparsers.add_parser(
        "outliers_table",
        description="Takes a table of values and converts to a "
                    "table of outlier counts.",
    )
    outliers_table.add_argument(
        "values",
        type=_is_valid_file,
        help="File path to input values. Columns must be samples, "
             "genes must be sites or genes. Only .tsv and .csv accepted.",
    )
    outliers_table.add_argument(
        "--output_prefix",
        type=_check_output_prefix,
        default="outliers",
        help="Output prefix for writing files. Default outliers. ",
    )
    outliers_table.add_argument(
        "--iqrs",
        type=_check_positive,
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
        type=_is_valid_file,
        help="Annotation table with samples as rows and "
             "annotation labels as columns. ",
    )
    binarize.add_argument(
        "--output_prefix",
        type=_check_output_prefix,
        default="annotations",
        help="Output prefix for writing files. Default annotations. Suffix will be '.binarized.tsv'",
    )

    compare_groups = subparsers.add_parser(
        "compare_groups",
        description="Takes an annotation table and outlier count "
                    "table (output of outliers_table) and outputs qvalues from "
                    "a statistical test that looks for enrichment of outlier values in each "
                    "group in the annotation table. For each value in each comparison, the qvalue table will "
                    "have 1 column, if there are any genes in that comparison. "
    )
    compare_groups.add_argument(
        "outliers_table",
        type=_is_valid_file,
        help="Table of outlier counts (output of outliers_table). Must be "
             ".tsv or .csv file, with outlier and non-outlier counts as columns, and genes/sites as "
             "rows. ",
    )
    compare_groups.add_argument(
        "annotations",
        type=_is_valid_file,
        help="Table of annotations. Must be .csv or .tsv. Samples as rows "
             "and comparisons as columns. Comparisons must have only  "
             "unique values (not including missing values). If there are "
             "more options than that, you can use binarize to prepare the "
             "table. ",
    )
    compare_groups.add_argument(
        "--ind_subset",
        type=_is_valid_file,
        default=None,
        help="File with subset of indexes to consider in comparison",
    )
    compare_groups.add_argument(
        "--ind_sep",
        type=str,
        default=None,
        help="Index separator for subsetting genes. Only needed if using ind_subset, and if rows "
             "of outliers are NOT aggregated. ",
    )
    compare_groups.add_argument(
        "--output_prefix",
        type=_check_output_prefix,
        default="outliers",
        help="Output prefix for writing files. Default outliers. ",
    )
    compare_groups.add_argument(
        "--frac_filter",
        type=_bn0and1,
        default=0.3,
        help="The minimum fraction of samples per group that must have an outlier in a gene to"
             "consider that gene in the analysis. This is used to prevent a high number of outlier "
             "values in 1 sample from driving a low qvalue. Default 0.3",
    )
    compare_groups.add_argument(
        "--write_comparison_summaries",
        default=False,
        action="store_true",
        help="Use flag to write a separate file for each column in the annotations table, "
             "with outlier counts in each group, p-values and q-values in each group. ",
    )
    compare_groups.add_argument(
        "--iqrs",
        type=_check_positive,
        default=None,
        help="Number of IQRs used to define outliers in the input count table. Optional.",
    )
    compare_groups.add_argument(
        "--up_or_down",
        type=str,
        choices=["up", "down"],
        help="Whether input outlier table represents up or down outliers. Needed for "
             "output file labels. Default up",
    )
    compare_groups.add_argument(
        "--write_gene_list",
        default=False,
        action="store_true",
        help="Use flag to write a list of significantly enriched genes for each value in each "
             "comparison. If used, need an fdr threshold as well. ",
    )
    compare_groups.add_argument(
        "--make_heatmaps",
        default=False,
        action="store_true",
        help="Use flag to draw a heatmap of signficantly enriched genes for each value in each "
             "comparison. If used, need an fdr threshold as well. ",
    )
    compare_groups.add_argument(
        "--fdr",
        type=_bn0and1,
        default=0.05,
        help="FDR cut off to use for signficantly enriched gene lists and heatmaps. Default 0.05",
    )
    compare_groups.add_argument(
        "--red_or_blue",
        type=str,
        choices=["red", "blue"],
        default="red",
        help="If --make_heatmaps is called, color of values to draw on heatmap. Default red. ",
    )
    compare_groups.add_argument(
        "--annotation_colors",
        type=_is_valid_file,
        default=None,
        help="File with color map to use for annotation header if --make_heatmaps is used. Must "
             "have a 'value    color' format for each value in annotations. Any value not "
             "represented will be assigned a new color. ",
    )

    visualize = subparsers.add_parser(
        "visualize",
        description="Used to make custom heatmaps from significant " "genes. ",
    )
    visualize.add_argument(
        "comparison_qvalues",
        type=_is_valid_file,
        help="Table of qvalues, output from compare_groups. Must be .csv or .tsv. Has genes/sites "
             "as rows and comparison values as columns. ",
    )
    visualize.add_argument(
        "annotations",
        type=_is_valid_file,
        help="Table of annotations used to generate qvalues. ",
    )
    visualize.add_argument(
        "visualization_table",
        type=_is_valid_file,
        help="Values to visualize in heatmap. Samples as columns and "
             "genes/sites as rows. Using outlier fraction table is recommended, but original "
             "values can also be used if no aggregation was used. ",
    )
    visualize.add_argument(
        "comparison_of_interest",
        type=str,
        help="Name of column in qvalues table from which to visualize significant genes. ",
    )
    visualize.add_argument(
        "--output_prefix",
        type=_check_output_prefix,
        default="outliers",
        help="Output prefix for writing files. Default outliers. ",
    )
    visualize.add_argument(
        "--annotations_to_show",
        type=str,
        default=None,
        nargs="+",
        help="Names of columns from the annotation table to show in the header of the heatmap. "
             "Default is all columns. ",
    )
    visualize.add_argument(
        "--fdr",
        type=_bn0and1,
        default=0.05,
        help="FDR threshold to use to select genes to visualize. Default 0.05",
    )
    visualize.add_argument(
        "--red_or_blue",
        type=str,
        choices=["red", "blue"],
        default="red",
        help="Color of values to draw on heatmap. Default red. ",
    )
    visualize.add_argument(
        "--annotation_colors",
        type=_is_valid_file,
        default=None,
        help="File with color map to use for annotation header. Must "
             "have a line with 'value    color' format for each value in annotations. Any value "
             "not represented will be assigned a new color. ",
    )
    visualize.add_argument(
        "--write_gene_list",
        default=False,
        action="store_true",
        help="Use flag to write a list of significantly enriched genes for each value in each "
             "comparison.",
    )

    deva = subparsers.add_parser(
        "deva",
        description="Runs whole outliers pipeline. Has options to output every "
                    "possible output. ",
    )
    deva.add_argument(
        "values",
        type=_is_valid_file,
        help="File path to input values. Samples are columns and genes/sites are rows. Only .tsv "
             "and .csv accepted.",
    )
    deva.add_argument(
        "annotations",
        type=_is_valid_file,
        help="File path to annotation values. Rows are sample names, "
             "header is different annotations. e.g. mutation status.",
    )
    deva.add_argument(
        "--output_prefix",
        type=_check_output_prefix,
        default="outliers",
        help="Output prefix for writing files. Default outliers. ",
    )
    deva.add_argument(
        "--iqrs",
        type=_check_positive,
        default=1.5,
        help="Number of inter-quartile ranges (IQRs) above or below the "
             "median to consider a value an outlier. Default is 1.5.",
    )
    deva.add_argument(
        "--up_or_down",
        type=str,
        default="true",
        choices=["up", "down"],
        help="Whether to look for up or down outliers. Choices are up or down. Default up.",
    )
    deva.add_argument(
        "--do_not_aggregate",
        default=False,
        action="store_true",
        help="Use flag if you do not want to sum outliers based on site prefixes.",
    )
    deva.add_argument(
        "--write_outlier_table",
        default=False,
        action="store_true",
        help="Use flag to write a table of outlier counts.",
    )
    deva.add_argument(
        "--write_frac_table",
        default=False,
        action="store_true",
        help="Use flag if you want to write a table with fraction of "
             "values per site per sample that are outliers. Useful for custom visualization.",
    )
    deva.add_argument(
        "--ind_sep",
        type=str,
        default="-",
        help="If site labels have a parent molecule (e.g. a gene name such "
             "as ATM) and a site identifier (e.g. S365) this is the "
             "delimiter between the two elements. Default is -",
    )
    deva.add_argument(
        "--frac_filter",
        type=_bn0and1,
        default=0.3,
        help="The minimum fraction of samples per group that must have an outlier in a gene to"
             "consider that gene in the analysis. This is used to prevent a high number of outlier "
             "values in 1 sample from driving a low qvalue. Default 0.3",
    )
    deva.add_argument(
        "--write_comparison_summaries",
        default=False,
        action="store_true",
        help="Use flag to write a separate file for each column in the annotations table, "
             "with outlier counts in each group, p-values and q-values in each group. ",
    )
    deva.add_argument(
        "--fdr",
        type=_bn0and1,
        default=0.05,
        help="FDR threshold to use to select genes to visualize. Default 0.05",
    )
    deva.add_argument(
        "--write_gene_list",
        default=False,
        action="store_true",
        help="Use flag to write a list of significantly enriched genes for each value in each "
             "comparison.",
    )
    deva.add_argument(
        "--make_heatmaps",
        default=False,
        action="store_true",
        help="Use flag to draw a heatmap of significantly enriched genes for each value in each "
             "comparison. If used, need an fdr threshold as well. ",
    )
    deva.add_argument(
        "--red_or_blue",
        type=str,
        choices=["red", "blue"],
        default="red",
        help="Color of values to draw on heatmap. Default red. ",
    )
    deva.add_argument(
        "--annotation_colors",
        type=_is_valid_file,
        default=None,
        help="File with color map to use for annotation header. Must "
             "have a line with 'value    color' format for each value in annotations. Any value "
             "not represented will be assigned a new color. ",
    )

    simulations = subparsers.add_parser(
        "simulations",
        description="Add here. ",
    )
    simulations.add_argument(
        "values",
        type=_is_valid_file,
        help="File path to input values. Samples are columns and genes/sites are rows. Only .tsv "
             "and .csv accepted.",
    )
    simulations.add_argument(
        "--ind_sep",
        type=str,
        default="-",
        help="Delimiter between the parent molecule (e.g. a gene name such "
             "as ATM) and a site identifier (e.g. S365). Default is -",
    )
    simulations.add_argument(
        "--iqrs",
        type=_check_positive,
        default=1.5,
        help="Number of inter-quartile ranges (IQRs) above or below the "
             "median to consider a value an outlier. Default is 1.5.",
    )
    simulations.add_argument(
        "--reps",
        type=_check_positive,
        default=1000000,
        help="Number of repetitions for the simulation to perform. "
             "Default is 1,000,000.",
    )
    simulations.add_argument(
        "--output_prefix",
        type=_check_output_prefix,
        default="simulated_pvals",
        help="Output prefix for writing files. Default is 'simulated_pvals'.",
    )
    simulations.add_argument(
        "--molecules",
        nargs='+',
        default=[],
        help="List of parent molecules of interest. Empty list or absence of "
             "argument defaults to all parent molecules in input file.",
    )
    simulations.add_argument(
        "--pval",
        type=_bn0and1,
        default=0.05,
        help="p-value threshold for significant results. Must be between 0 and 1."
        "Default is 0.05.",
    )

    return parser


# Module runner
def _main(args: Optional[List[str]] = None):
    if args is None:
        args = sys.argv[1:]
    args = _make_parser().parse_args(args)

    logger = _set_up_logger(args.output_prefix)

    logger.info("Running deva in %s mode" % args.which)
    for arg in vars(args):
        logger.info("Parameter %s: %s" % (arg, getattr(args, arg)))

    if args.which == "outliers_table":
        df = parsers.read_in_values(args.values)
        make_outliers_table(
            df,
            iqrs=args.iqrs,
            up_or_down=args.up_or_down,
            aggregate=~args.do_not_aggregate,
            save_outlier_table=True,
            save_frac_table=args.write_frac_table,
            output_prefix=args.output_prefix,
            ind_sep=args.ind_sep,
        )

    elif args.which == "binarize":
        annotations = parsers.read_in_values(args.annotations)
        annotations = parsers.binarize_annotations(annotations)
        annotations.to_csv("%s.binarized.tsv" % args.output_prefix, sep="\t")

    elif args.which == "normalize":
        target = parsers.read_in_values(args.unnormed_values)
        df = parsers.normalize(target)
        df.to_csv(args.output_prefix + ".normalized.tsv", sep='\t')

    elif args.which == "compare_groups":
        outliers = parsers.read_in_outliers(
            args.outliers_table, args.up_or_down, args.iqrs
        )
        if args.ind_subset:
            with open(args.ind_subset, 'r') as fh:
                ind_list = [i.strip() for i in fh.readlines()]
                outliers.df = subset_by_genes(outliers.df, ind_list, args.ind_sep)

        annotations = parsers.read_in_values(args.annotations)
        qVals = compare_groups_outliers(
            outliers,
            annotations,
            frac_filter=args.frac_filter,
            save_qvalues=True,
            output_prefix=args.output_prefix,
            save_comparison_summaries=args.write_comparison_summaries,
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
        qvals = parsers.read_in_values(args.comparison_qvalues)
        annotations = parsers.read_in_values(args.annotations)
        frac_table = parsers.read_in_values(args.visualization_table)
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

    elif args.which == "deva":
        df = parsers.read_in_values(args.values)
        annotations = parsers.read_in_values(args.annotations)
        Outliers, qVals = deva(
            df,
            annotations,
            iqrs=args.iqrs,
            frac_filter=args.frac_filter,
            up_or_down=args.up_or_down,
            aggregate=~args.do_not_aggregate,
            save_outlier_table=args.write_outlier_table,
            save_frac_table=args.write_frac_table,
            save_qvalues=True,
            output_prefix=args.output_prefix,
            ind_sep=args.ind_sep,
            save_comparison_summaries=args.write_comparison_summaries,
        )
        if args.write_gene_list:
            qVals.write_gene_lists(args.fdr, args.output_prefix)

        if args.make_heatmaps:
            for col_of_interest in qVals.df.columns:
                plot_heatmap(
                    annotations,
                    qVals.df,
                    col_of_interest,
                    Outliers.frac_table,
                    fdr=args.fdr,
                    red_or_blue=args.red_or_blue,
                    output_prefix=args.output_prefix,
                    colors=args.annotation_colors,
                    savefig=True,
                )
                plt.close()

    elif args.which == "simulations":
        run_simulations(
            args.values,
            args.ind_sep,
            args.iqrs,
            int(args.reps),
            args.output_prefix,
            args.molecules,
            args.pval,
        )

    with open(parameters_file_name % args.output_prefix, "w") as fh:
        for arg in vars(args):
            fh.write("%s: %s\n" % (arg, getattr(args, arg)))


# Run cli
if __name__ == "__main__":
    _main()
