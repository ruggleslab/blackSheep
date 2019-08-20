from typing import Optional, Iterable
import logging
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import seaborn as sns
from blacksheep import catheat
from blacksheep._constants import *


def _get_sample_order(annotations: DataFrame, col_of_interest: str) -> Iterable:
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


def _get_genes(qvals: DataFrame, fdr: float, col: str) -> list:
    """Collects significant genes

    Args:
        qvals: qvalues DataFrame
        fdr: FDR cut off
        col: Column for which to collect genes

    Returns:
        List of significant genes

    """

    return list(qvals.loc[(qvals[col] < fdr), :].index)


def _pick_color(red_or_blue: str):
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


def _check_colors(colors: dict) -> dict:
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


def _assign_colors(data: DataFrame, cmap: dict, palette) -> dict:
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
    colors = catheat._gen_colors(palette, len(missing_entries))
    cmap.update({v: colors[i] for i, v in enumerate(missing_entries)})

    return cmap


def _determine_colors(path: str, annotations: DataFrame) -> dict:
    """Takes a file with value, color pairs and fills out any other needed colors.

    Args:
        path: File path to value, color pairs
        annotations: Annotation DataFrame

    Returns:
        Color dictionary

    """

    if not path:
        return _assign_colors(annotations, {}, default_palette)
    try:
        with open(path, "r") as fh:
            colors = {line.split()[0]: line.split()[1] for line in fh.readlines()}
        colors = _check_colors(colors)
    except FileNotFoundError:
        logging.warning("%s is not a valid file, generating colors" % path)
        colors = {}

    return _assign_colors(annotations, colors, default_palette)


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
    sample_order = _get_sample_order(annotations, annot_label)
    genes = _get_genes(qvals, fdr, col_of_interest)
    if not genes:
        logging.warning("No significant genes at FDR %s in %s" % (fdr, col_of_interest))
        return None

    annotations = annotations.reindex(sample_order, axis=1)
    vis_table = vis_table.reindex(genes).reindex(sample_order, axis=1)

    # Get colors
    cmap = _pick_color(red_or_blue)
    colors = _determine_colors(colors, annotations)

    # Get label
    label = col_of_interest[10:]

    # Set up figure
    sns.set(font="arial", style="white", color_codes=True, font_scale=1)
    plot_height = min(max((0.19 * (len(annotations) + len(genes))), 2), 11)
    plot_width = min(max((0.15 * len(annotations.columns)), 4), 8.5)
    fig = plt.figure(figsize=(plot_width, plot_height))
    gs = plt.GridSpec(
        figure=fig,
        nrows=3,
        ncols=2,
        width_ratios=[len(annotations.columns), 2],
        height_ratios=[len(annotations)] + [len(vis_table) / 2 for i in range(0, 2)],
        wspace=0.01,
        hspace=0.01,
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
        leg_ax=leg_ax,
        leg_kws=dict(loc=(0, 0.05), facecolor="white", edgecolor="white"),
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

    if savefig:
        fig_path = figure_file_name % (output_prefix, label, fdr)
        logging.info("Saving figure to %s" % fig_path)
        plt.savefig(fig_path, dpi=200, bbox_inches="tight")

    return [annot_ax, vals_ax, cbar_ax, leg_ax]
