#!/usr/bin/env python
# coding: utf-8



import anndata as ad
import seaborn as sns
import pandas as pd
import numpy as np
from pycdr.utils import get_top_genes
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects


# specific details for scatterplot
#['factor.1002', 'factor.1553', 'factor.976']
coords_1002 = [(0,20), (0,20), (0,20), (0,20), (0,20)]
coords_1553 = [(0,20), (40,0), (25,25), (-25,20), (0,-30)]
coords_976 = [(0,20), (40,0), (40,0), (-25,20), (0,-30)]


def process_data_for_barplots(x, factor_of_interest, degenes, geneindex):
    
    genesofinterest = get_top_genes(x, factor_of_interest).sort_values("z_score", ascending = False)
    genes = genesofinterest.head(20).index.to_list()
    
    loadingscores = get_top_genes(x, factor_of_interest).sort_values("z_score", ascending = False)
    scores = loadingscores.head(20).z_score.to_list()
    color_list = []
    
    genes_plotting = pd.DataFrame(zip(genes, scores), index = genes).reindex(geneindex)

    Product = genes_plotting[0]

    for i in Product:
        if i in degenes.to_list():
            color_list.append("black")
        else:
            color_list.append("white")

    Quantity = genes_plotting[1]
    
    return Product, Quantity, color_list


def plot_barplots(Product, Quantity, color_list, filename):

    fig, ax = plt.subplots(figsize=(2,5))

    # Horizontal Bar Plot


    ax.barh(Product, Quantity, color = color_list, edgecolor = "black")
    ax.invert_yaxis()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)

    import matplotlib.ticker as ticker
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.grid(axis='x')

    #plt.ylabel('Gene')
    plt.xlabel('Z-scores', fontsize = "x-large")
    plt.tight_layout()
    plt.savefig(filename, dpi = 300)



def process_data_for_scatterplot(x, factor, genelist, pheno1, pheno2):
    
    ctrl = x[x.obs[factor] == pheno1].copy()
    stim = x[x.obs[factor] == pheno2].copy()
    genesstim = sc.pp.calculate_qc_metrics(stim)[1]
    genesctrl = sc.pp.calculate_qc_metrics(ctrl)[1]
    factor0 = genelist
    
    stimcounts = genesstim[genesstim.index.isin(factor0)].log1p_mean_counts
    nonstimcounts = genesctrl[genesctrl.index.isin(factor0)].log1p_mean_counts

    x1 = stimcounts
    y1 = nonstimcounts
    return stimcounts, nonstimcounts



def process_data_for_heatmaps(x, factor_of_interest):
    genes = get_top_genes(x, factor_of_interest).sort_values("z_score", ascending = False)
    genesofinterest = genes.head(20).index.to_list() 

    inder = x.var.features.isin(genesofinterest)

    responder_adata = x[x.obs.stim == "STIM",inder]
    nonresponder_adata = x[x.obs.stim == "CTRL",inder]

    responder_adata_mat = responder_adata.X.toarray()
    nonresponder_adata_mat = nonresponder_adata.X.toarray()

    corr = np.corrcoef(responder_adata_mat.T)
    corr2 = np.corrcoef(nonresponder_adata_mat.T)

    resp = pd.DataFrame(corr)
    resp.columns = responder_adata.var.index
    resp.index = responder_adata.var.index

    nonresp = pd.DataFrame(corr2)
    nonresp.columns = nonresponder_adata.var.index
    nonresp.index = nonresponder_adata.var.index
    
    p = sns.clustermap(resp)

    order = p.dendrogram_row.reordered_ind

    geneindex = resp.columns[order]

    nonresp = nonresp.reindex(geneindex, axis=1)
    nonresp = nonresp.reindex(geneindex, axis=0)

    resp = resp.reindex(geneindex, axis=1)
    resp = resp.reindex(geneindex, axis=0)
    
    return resp, nonresp, genesofinterest



def plot_heatmaps(resp, nonresp, filename, thres = 0.10):
    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (10,5), gridspec_kw={'width_ratios': [1, 1]})

    sns.heatmap(
        resp, 
        square=True,
        cmap = 'RdBu_r',
        vmax = thres,
        ax = ax1,
        cbar = None,
    )

    sns.heatmap(
        nonresp, 
        square= True,
        cmap = 'RdBu_r',
        vmax = thres,
        ax = ax2,
        cbar = None,
        yticklabels=False

    )

    ax1.set_title("IFN-stimulated", fontsize = "x-large")
    ax2.set_title("Control", fontsize = "x-large")

    plt.tight_layout()

    plt.savefig(filename, dpi = 300)



def plot_scatterplot(X, Y, filename):

    fig, ax = plt.subplots(figsize = (5,5))
    #ax.plot(x, y, ls='--', color = 'red')

    X = a
    Y = b
    #X = a.reindex(genes)
    #Y = b.reindex(genes)

    toplevel = np.ceil(np.max([X,Y]))

    toplevel = np.round(np.max([X,Y]), 1)

    plt.ylim(-0.05, toplevel + 0.05)
    plt.xlim(-0.05, toplevel + 0.05)

    plt.xticks(np.arange(0, toplevel + 0.05, 0.1))
    plt.yticks(np.arange(0, toplevel + 0.05, 0.1))

    ax.grid()

    I = 5

    Px = X[0:I]
    Py = Y[0:I]

    ax.scatter(X, Y, edgecolor="None", facecolor="gray", alpha=0.5)
    ax.scatter(Px, Py, edgecolor="black", facecolor="white", zorder=20)
    ax.scatter(Px, Py, edgecolor="black", facecolor="C1", alpha=0.5, zorder=30)

    plt.plot([0, toplevel], [0, toplevel], ls='--', color = 'red') # plots line y = x

    y, dy = 1.0, 0.125
    style = "arc,angleA=-0,angleB=0,armA=-100,armB=0,rad=0"

    GENENAMES = X.index.to_list()

    for i in range(I):
        text = ax.annotate(
            GENENAMES[i],
            xy=(Px[i], Py[i]),
            xycoords="data",
            xytext=(0, 20),
            textcoords="offset points",
            ha="center",
            size="large",
            arrowprops=dict(
                arrowstyle="->", shrinkA=0, shrinkB=5, color="black", linewidth=0.75
            ),
        )
        text.set_path_effects(
            [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
        )
        text.arrow_patch.set_path_effects(
            [path_effects.Stroke(linewidth=2, foreground="white"), path_effects.Normal()]
        )

    plt.xlabel("log1p expression: CTRL monocytes", fontsize = "large")
    plt.ylabel("log1p expression: STIM monocytes", fontsize = "large")

    plt.title("Key genes in factor loading", fontsize = "x-large")

    #ax.spines['top'].set_linestyle((0, (10, 10)))
    #ax.spines['right'].set_linestyle((0, (10, 10, 1, 10)))
    plt.savefig(filename, dpi = 300)



x = ad.read(snakemake.input[0])

degenes = pd.read_csv(snakemake.input[1])["0"]


a,b,c = process_data_for_heatmaps(x, 976)
plot_heatmaps(a,b, "results/heatmaps_976.png", 0.15)
e,f,g = process_data_for_barplots(x, 976, degenes, a.index)
plot_barplots(e,f,g, "results/barplot_976.png")


a,b,c = process_data_for_heatmaps(x, 1002)
plot_heatmaps(a,b, "results/heatmaps_1002.png", 0.15)
e,f,g = process_data_for_barplots(x, 1002, degenes, a.index)
plot_barplots(e,f,g, "results/barplot_1002.png")


a,b,c = process_data_for_heatmaps(x, 1553)
plot_heatmaps(a,b, "results/heatmaps_1553.png", 0.15)
e,f,g = process_data_for_barplots(x, 1553, degenes, a.index)
plot_barplots(e,f,g, "results/barplot_1553.png")

