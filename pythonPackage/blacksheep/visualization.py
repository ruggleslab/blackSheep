from typing import Optional
from pandas import DataFrame
import numpy as np
import catheat
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import seaborn as sns


def getSampleOrder(annotations: DataFrame, col_of_interest: str):
    sort_by = [col for col in annotations.index if col != col_of_interest]
    annotations = annotations.sort_values([col_of_interest] + sort_by, axis=1)
    return annotations.columns


def getGenes(qvals: DataFrame, fdr: float, col: str) -> list:
    return list(qvals.loc[(qvals[col] < fdr), :].index)


def pickColor(red_or_blue: str):
    if red_or_blue == "red":
        cmap = sns.cubehelix_palette(
            start=0.857,
            rot=0.00,
            gamma=1.5,
            hue=1,
            light=1,
            dark=0.2,
            reverse=False,
            as_cmap=True,
        )
    elif red_or_blue == "blue":
        cmap = sns.cubehelix_palette(
            start=3,
            rot=0.00,
            gamma=1.5,
            hue=1,
            light=1,
            dark=0.2,
            reverse=False,
            as_cmap=True,
        )
    else:
        # TODO log error here
        raise ValueError("Invalid color choice, must be red or blue, setting to red.")

    cmap.set_bad("#BDBDBD")
    return cmap


def checkColors(colors: dict) -> dict:
    for lab, color in colors.items():
        try:
            mpatch.Patch(color=color)
        except ValueError:
            # TODO change to log
            print("%s is not a valid color" % color)
            colors.pop(lab)
    return colors


def _gen_colors(pal, n: int):
    """ Generate colours from provided palette.
    """

    # If string
    if isinstance(pal, str):
        # First, check if seaborn palette
        try:
            colors = sns.color_palette(pal, n)
        # If not, try getting the matplotlib palette
        except:
            pal = plt.get_cmap(pal)
            colors = [pal(i) for i in np.linspace(0, 1, n)]
    # If palette provided
    elif isinstance(pal, mcolors.LinearSegmentedColormap):
        colors = [pal(i) for i in np.linspace(0, 1, n)]
    # If list of colors
    elif isinstance(pal, (list, np.ndarray)):
        if len(pal) < n:
            raise ValueError(
                "Must provide at least as many colors as there are unique entries: {0}".format(
                    len(unique_values)
                )
            )
        else:
            colors = pal
    else:
        raise TypeError(
            'Unable to generate colors from palette of type "{0}"'.format(type(pal))
        )

    return colors


def listToFile(lis: list, filename: str):
    with open(filename, "w") as fh:
        for x in lis:
            fh.write("%s\n" % x)


def assign_colors(data: DataFrame, cmap: dict, palette) -> dict:

    unique_values = sorted(np.unique(data.values.astype(str)))
    n_unique = len(unique_values)

    missing_entries = [v for v in unique_values if v not in cmap.keys()]
    colors = _gen_colors(palette, len(missing_entries))

    cmap.update({v: colors[i] for i, v in enumerate(missing_entries)})
    return cmap


def parseColors(path: str, annotations: DataFrame) -> dict:
    if path is None:
        colors = {}
    else:
        try:
            with open(path, "r") as fh:
                colors = {line.split()[0]: line.split()[1] for line in fh.readlines()}
            colors = checkColors(colors)

        except FileNotFoundError:
            # TODO log
            print("%s is not a valid file" % path)
            colors = {}

    colors = assign_colors(annotations, colors, "Set2")
    return colors


def plot_heatmap(
    annotations: DataFrame,
    qvals: DataFrame,
    col_of_interest: str,
    vis_table: DataFrame,
    fdr: float = 0.05,
    red_or_blue: str = "red",
    output_prefix: str = "outliers",
    colors: Optional[str] = None,
    savefig: bool = True,
) -> list:
    """
    TODO fill out help
    :param annotations:
    :param qvals:
    :param col_of_interest:
    :param vis_table:
    :param fdr:
    :param red_or_blue:
    :param output_prefix:
    :param colors:
    :param savefig:
    :return:
    """

    annotations = annotations.transpose()
    vis_table = vis_table.transpose()

    # Get orders
    annot_label = col_of_interest.split("_", 1)[1].rsplit("_", 1)[0]
    sample_order = getSampleOrder(annotations, annot_label)
    genes = getGenes(qvals, fdr, col_of_interest)
    if len(genes) == 0:
        # TODO logging, add a error here
        print("No signficant genes at %s" % fdr)
        return None

    annotations = annotations.reindex(sample_order, axis=1)
    vis_table = vis_table.reindex(genes).reindex(sample_order, axis=1)

    # Get colors
    cmap = pickColor(red_or_blue)
    colors = parseColors(colors, annotations)

    # Get label
    label = col_of_interest[10:]

    # Set up figure
    sns.set(font="arial", style="white", color_codes=True, font_scale=0.60)
    plot_height = min(max(0.1 * (len(annotations) + len(genes)), 1), 15)
    plot_width = 0.05 * len(annotations.columns)

    fig = plt.figure(figsize=(plot_width, plot_height))
    gs = plt.GridSpec(
        figure=fig,
        nrows=3,
        ncols=2,
        width_ratios=[len(annotations.columns), 2],
        height_ratios=[len(annotations)] + [len(vis_table) / 2 for i in range(0, 2)],
        wspace=0.01,
        hspace=0.01,
        left=0.20,
        right=0.80,
        top=0.88,
        bottom=0.05,
    )

    annot_ax = plt.subplot(gs[0, 0])
    vals_ax = plt.subplot(gs[1:, 0])
    cbar_ax = plt.subplot(gs[-1, 1])
    leg_ax = plt.subplot(gs[:-1, 1])
    leg_ax.axis("off")

    # Header
    catheat.heatmap(
        annotations,
        cmap=colors,
        ax=annot_ax,
        legend=False,
        xticklabels=False,
        yticklabels=annotations.index,
    )

    annot_ax.set_title("Outliers in %s" % col_of_interest, va="top")
    annot_ax.set_yticklabels(annotations.index, rotation=0)

    # Values
    sns.heatmap(
        vis_table,
        ax=vals_ax,
        cbar_ax=cbar_ax,
        cmap=cmap,
        vmin=0,
        vmax=1,
        yticklabels=vis_table.index,
        xticklabels=False,
        cbar_kws=dict(label="Fraction\nOutliers"),
    )
    vals_ax.set_xlabel("")
    vals_ax.set_ylabel("")
    vals_ax.set_yticklabels(vis_table.index, rotation=0)

    # Legend
    handles = [mpatch.Patch(label=l, color=c) for l, c in colors.items()]
    leg_ax.legend(
        handles=handles,
        loc=(0, 0.05),
        facecolor="white",
        edgecolor="white",
        labelspacing=0,
    )
    if savefig:
        plt.savefig("%s.%s.heatmap.pdf" % (output_prefix, label), dpi=200)

    return [annot_ax, vals_ax, cbar_ax, leg_ax]


def write_genes(qvals: DataFrame, fdr: float, output_prefix: str):
    for col in qvals.columns:
        genes = getGenes(qvals, fdr, col)
        if genes:
            listToFile(genes, "%s.fdr%s.%s.txt" % (output_prefix, fdr, col))
