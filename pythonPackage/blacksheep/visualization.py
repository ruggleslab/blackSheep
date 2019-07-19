from typing import Optional, List
from pandas import DataFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


def computeOrder(df, optimal=True, dist_method="euclidean", cluster_method="average"):

    dist_mat = pdist(df, metric=dist_method)
    link_mat = hierarchy.linkage(dist_mat, method=cluster_method)

    if optimal == True:
        return hierarchy.leaves_list(
            hierarchy.optimal_leaf_ordering(link_mat, dist_mat)
        )
    else:
        return hierarchy.leaves_list(link_mat)


def clustermap(
    df,
    dist_method="euclidean",
    cluster_method="average",
    col_cluster=True,
    row_cluster=True,
    optimal=True,
    **heatmap_kws
):

    if row_cluster == True:
        row_order = computeOrder(df, optimal, dist_method, cluster_method)
        row_order = [df.index[i] for i in row_order]
    else:
        row_order = df.index

    if col_cluster == True:
        col_order = computeOrder(df.transpose(), optimal, dist_method, cluster_method)
        col_order = [df.columns[i] for i in col_order]
    else:
        col_order = df.columns

    df = df.reindex(col_order, axis=1).reindex(row_order, axis=0)

    ax = sns.heatmap(df, **heatmap_kws)

    return ax, df


def plot_heatmap():
    pass
