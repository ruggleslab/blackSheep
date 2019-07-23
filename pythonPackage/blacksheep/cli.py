import sys
import argparse
import os.path
from .outliers import run_outliers
from .outliers import make_outliers_table
from .outliers import compare_groups_outliers
from .classes import OutlierTable
from . import parsers


def check_positive(arg):
    try:
        arg = float(arg)
    except:
        raise argparse.ArgumentTypeError("%s is an invalid positive value" % arg)
    if arg < 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive value" % arg)
    else:
        return arg


def is_valid_file(arg):
    if not os.path.exists(arg):
        raise argparse.ArgumentTypeError("The file %s does not exist" % arg)
    elif (arg[-4:] != ".tsv") and (arg[-4:] != ".csv"):
        raise argparse.ArgumentTypeError("%s is not a tsv or csv" % arg)
    else:
        return arg


def check_output_prefix(arg):
    if "/" in arg:
        prefix = arg.rsplit("/", 1)[0]
        if not os.path.exists(prefix):
            raise argparse.ArgumentTypeError("%s does not exist" % prefix)
        else:
            return arg
    else:
        return arg


def bn0and1(arg):
    try:
        arg = float(arg)
    except:
        raise argparse.ArgumentTypeError("%s is an invalid positive value" % arg)
    if (type(arg) != float) or (arg < 0) or (arg > 1):
        raise argparse.ArgumentTypeError("%s is not a number between 0 and 1" % arg)
    else:
        return arg

def parse_args(args):
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

    return parser.parse_args(args)


def main(args):
    if args is None:
        args = parse_args(sys.argv[1:])
    else:
        args = parse_args(args)
    if args.which == "outliers_table":
        df = parsers.parseValues(args.values)
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

    elif args.which == "compare_groups":
        outliers = parsers.parseOutliers(
            args.outliers_table, args.up_or_down, args.iqrs, None, None
        )
        outliers.samples = [ind.rsplit("_", 1)[0] for ind in outliers.df.index]

        annotations = parsers.parseAnnotations(args.annotations_table)
        qvalues = compare_groups_outliers(
            outliers,
            annotations,
            args.frac_filter,
            save_qvalues=True,
            output_prefix=args.output_prefix,
            output_comparison_summaries=args.write_comparison_summaries,
        )

    elif args.which == "outliers":
        df = parsers.parseValues(args.values)
        annotations = parsers.parseAnnotations(args.annotations)
        qvals, frac_table = run_outliers(
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
    with open("%s.parameters.txt" % args.output_prefix, "w") as fh:
        for arg in vars(args):
            fh.write("%s: %s\n" % (arg, getattr(args, arg)))


if __name__ == "__main__":
    args = None
    main(args)
