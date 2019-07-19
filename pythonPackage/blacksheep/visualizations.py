from pandas import DataFrame
import pandas as pd
import numpy as np

import catheat
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


def clustermap(
        df: DataFrame,
        fillna: float = 0,
        dist_method: str = "euclidean",
        cluster_method: str = "average",
        col_cluster: bool = True,
        row_cluster: bool = True,
        optimal: bool = True,
        **heatmap_kws
):
    df = df.fillna(fillna).transpose()
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
    
    return
