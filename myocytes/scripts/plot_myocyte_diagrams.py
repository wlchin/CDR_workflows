import pandas as pd
import anndata as ad
import scanpy as sc
from pycdr.utils import get_top_genes
from pycdr.perm import calculate_enrichment
import matplotlib.ticker as ticker
import numpy as np
import matplotlib.pyplot as plt

adata = ad.read(snakemake.input[0])

factors = ["factor.1", "factor.5", "factor.54","factor.0"]
factor_nums = [1,5,54,0]

fig, ax = plt.subplots(1,4, figsize = (9,3))
for i in range(4):
    factor_loading = factors[i]
    info = adata.uns["dict_res_prop"][factor_loading]
    timepoints = np.array(info[2]).astype("str")
    vals = (np.array(info[0]).astype("int"))/(np.array(info[1]).astype("int")) * 100
    ax[i].bar(timepoints,vals, color = "black")
    ax[i].set_xlabel("Hours", fontsize = "x-large")
    ax[i].set_title(factor_loading, fontsize = "x-large", pad=10)
    ax[i].set_ylim((0,35))
    ax[i].spines['top'].set_visible(False)
    ax[i].spines['right'].set_visible(False)
    ax[i].spines['bottom'].set_visible(True)
    
fig.supylabel("% activated cells", fontsize = "xx-large")
plt.savefig("results/muscle_factors.png", dpi = 600)


fig, ax = plt.subplots(1,4, figsize = (9,6))
for i in range(4):
    factor_id = factor_nums[i]
    top20 = get_top_genes(adata, factor_id).sort_values("z_score", ascending = False).head(20)
    genesofinterest = top20.index.to_list()
    loadingscores = top20.z_score.to_list()

    ax[i].barh(genesofinterest, loadingscores, edgecolor = "black")
    ax[i].invert_yaxis()
    ax[i].spines['top'].set_visible(False)
    ax[i].spines['right'].set_visible(False)
    ax[i].spines['bottom'].set_visible(True)
    ax[i].spines['left'].set_visible(False)

    ax[i].xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax[i].grid(axis='x')
    ax[i].set_xlabel('Z-scores')

fig.supylabel("Top genes in factor loading", fontsize = "xx-large")
fig.tight_layout(pad = 4.0)
plt.savefig("results/muscle_genes.png", dpi = 300)


factors = ["factor.7", "factor.8"]

fig, ax = plt.subplots(1,2, figsize = (9,3))                                                                                                                                                                         
for i in range(2):                                                                                                                                                                                                   
    factor_loading = factors[i]                                                                                                                                                                                      
    info = adata.uns["dict_res_prop"][factor_loading]                                                                                                                                                                
    timepoints = np.array(info[2]).astype("str")                                                                                                                                                                     
    vals = (np.array(info[0]).astype("int"))/(np.array(info[1]).astype("int")) * 100                                                                                                                                 
    ax[i].bar(timepoints,vals, color = "black")                                                                                                                                                                      
    ax[i].set_xlabel("Hours", fontsize = "x-large")                                                                                                                                                                                                                                                                                                         
    ax[i].set_title(factor_loading, fontsize = "x-large", pad=10)                                                                                                                                                    
    ax[i].set_ylim((0,35))                                                                                                                                                                                           
    ax[i].spines['top'].set_visible(False)                                                                                                                                                                           
    ax[i].spines['right'].set_visible(False)                                                                                                                                                                         
    ax[i].spines['bottom'].set_visible(True)  
    
fig.supylabel("% activated cells", fontsize = "xx-large")
fig.savefig("results/bar_chart_subpopulations.png", dpi = 500)

