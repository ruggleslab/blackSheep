from typing import Optional, Iterable
import logging
from pandas import DataFrame
import numpy as np
import catheat
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import seaborn as sns
from .constants import *


def get_sample_order(annotations: DataFrame, col_of_interest: str) -> Iterable:
    """Orders the samples using annotations

    Args:
        annotations: comparisons vs samples annotation DataFrame
        col_of_interest: Comparison to sort by first

    Returns:
        Sorted order of samples

    """

    sort_by = [col for col in annotations.index if col != col_of_interest]
    annotations = annotations.sort_values([col_of_interest] + sort_by, axis=1)
    return annotations.columns


def get_genes(qvals: DataFrame, fdr: float, col: str) -> list:
    """Collects significant genes

    Args:
        qvals: qvalues DataFrame
        fdr: FDR cut off
        col: Column for which to collect genes

    Returns:
        List of significant genes

    """

    return list(qvals.loc[(qvals[col] < fdr), :].index)


def pick_color(red_or_blue: str):
    """Sets colormap for heatmap

    Args:
        red_or_blue: Use red or blue.

    Returns:
        Colormap

    """
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
        raise ValueError("Invalid color choice, must be red or blue, setting to red.")

    cmap.set_bad("#BDBDBD")
    return cmap


def check_colors(colors: dict) -> dict:
    """Makes sure every input color can be used as a color by the heatmap and legend.

    Args:
        colors: Dictionary of {values: colors}

    Returns:
        Dictionary of values: colors with invalid ones removed.

    """
    for lab, color in colors.items():
        try:
            mpatch.Patch(color=color)
        except ValueError:
            logging.warning("%s is not a valid color" % color)
            colors.pop(lab)
    return colors


def gen_colors(pal, n: int):
    """Generates colors from provided palette.

    Args:
        pal: List of colors, LinearSegmentedColormap or seaborn color palette. Must have
    enough colors for needed unique values.
        n: How many unique colors are needed

    Returns:
        List of colors to assign to values.

    """

    # If string
    if isinstance(pal, str):
        try:
            colors = sns.color_palette(pal, n)
        except:
            pal = plt.get_cmap(pal)
            colors = [pal(i) for i in np.linspace(0, 1, n)]
    elif isinstance(pal, mcolors.LinearSegmentedColormap):
        colors = [pal(i) for i in np.linspace(0, 1, n)]
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


def assign_colors(data: DataFrame, cmap: dict, palette) -> dict:
    """Combines provided colors, and adds more if needed

    Args:
        data: Annotations heatmap to color
        cmap: Provided colormap
        palette: Palette to use to generate unspecific colors

    Returns:
        Color dictionary for every unique value in annotations.

    """

    unique_values = sorted(np.unique(data.values.astype(str)))
    cmap = {v: c for v, c in cmap.items() if v in unique_values}
    missing_entries = [v for v in unique_values if v not in cmap.keys()]
    colors = gen_colors(palette, len(missing_entries))
    cmap.update({v: colors[i] for i, v in enumerate(missing_entries)})

    return cmap


def determine_colors(path: str, annotations: DataFrame) -> dict:
    """Takes a file with value, color pairs and fills out any other needed colors.

    Args:
        path: File path to value, color pairs
        annotations: Annotation DataFrame

    Returns:
        Color dictionary

    """

    if path is None:
        colors = {}
    else:
        try:
            with open(path, "r") as fh:
                colors = {line.split()[0]: line.split()[1] for line in fh.readlines()}
            colors = check_colors(colors)
        except FileNotFoundError:
            logging.warning("%s is not a valid file, generating colors" % path)
            colors = {}

    colors = assign_colors(annotations, colors, default_palette)
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
    savefig: bool = False,
) -> list:
    """Plots a heatmap of significantly enriched values for a given comparison.

    Args:
        annotations: Annotations DataFrame, samples as rows, annotations as columns
        qvals: qvalues DataFrame with genes/sites as rows and comparisons as columns
        col_of_interest: Which column from qvalues should be used to find signficant genes
        vis_table: Table to be visualized in heatmap. Index values should correspond to the \
        annotation df index, column names should correspond to qvals df index
        fdr: FDR threshold to for significance
        red_or_blue: Whether heatmap should be in red or blue color scale
        output_prefix:  If saving files, output prefix
        colors: File to find color map for annotation header
        savefig: Whether to save the plot to a pdf

    Returns: [annot_ax, vals_ax, cbar_ax, leg_ax]
        List of matplotlib axs, can be further customized before saving. In order the axes \
        contain: annotation header, the heatmap, the color bar, and the legend.

    """

    annotations = annotations.transpose()

    # Get orders
    annot_label = [col for col in annotations.index if col in col_of_interest][0]
    sample_order = get_sample_order(annotations, annot_label)
    genes = get_genes(qvals, fdr, col_of_interest)
    if not genes:
        logging.warning("No significant genes at FDR %s in %s" % (fdr, col_of_interest))
        return None

    annotations = annotations.reindex(sample_order, axis=1)
    vis_table = vis_table.reindex(genes).reindex(sample_order, axis=1)

    # Get colors
    cmap = pick_color(red_or_blue)
    colors = determine_colors(colors, annotations)

    # Get label
    label = col_of_interest[10:]

    # Set up figure
    sns.set(font="arial", style="white", color_codes=True, font_scale=1)
    plot_height = min(max((0.19 * (len(annotations) + len(genes))), 2), 11)
    plot_width = min(max((0.15 * len(annotations.columns)), 4), 8.5)
    margin_sizeLR = 1.75
    margin_sizeTB = 0.4
    fig = plt.figure(figsize=(plot_width, plot_height))
    gs = plt.GridSpec(
        figure=fig,
        nrows=3,
        ncols=2,
        width_ratios=[len(annotations.columns), 2],
        height_ratios=[len(annotations)] + [len(vis_table) / 2 for i in range(0, 2)],
        wspace=0.01,
        hspace=0.01,
        bottom=(margin_sizeTB / plot_height),
        top=1 - (margin_sizeTB / plot_height),
        left=0.05 + (margin_sizeLR / plot_width),
        right=1.05 - (margin_sizeLR / plot_width),
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

    annot_ax.set_title(plot_title % col_of_interest, fontsize=14)
    annot_ax.set_yticklabels(annotations.index, rotation=0)
    annot_ax.set_xlabel("")
    annot_ax.set_ylabel("")

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
        cbar_kws=dict(label=cbar_label),
    )
    vals_ax.set_yticklabels(vis_table.index, rotation=0)
    vals_ax.set_xlabel("")
    vals_ax.set_ylabel("")

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
        fig_path = figure_file_name % (output_prefix, label, fdr)
        logging.info("Saving figure to %s" % fig_path)
        plt.savefig(fig_path, dpi=200)

    return [annot_ax, vals_ax, cbar_ax, leg_ax]
