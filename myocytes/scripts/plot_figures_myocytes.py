#!/usr/bin/env python
# coding: utf-8

# In[81]:

import pandas as pd
import anndata as ad
import scanpy as sc
from pycdr.welford_gene_selection import get_df_for_factor
from pycdr.perm import get_df_loadings
from pycdr.perm import calculate_enrichment

# In[82]:


muscle = ad.read(snakemake.input[0])


# In[83]:


calculate_enrichment(muscle, "Hours", ["factor.1", "factor.0", "factor.8", "factor.5", "factor.7", "factor.54", "factor.8"], 100, "features", 0.1, 19)


# In[84]:




#sc.pp.highly_variable_genes(muscle, min_mean=0.0125, max_mean=3, min_disp=0.5)
#muscle = muscle[:, muscle.var.highly_variable]
sc.pp.scale(muscle, max_value=10)
sc.tl.pca(muscle, svd_solver='arpack')
sc.pp.neighbors(muscle, n_neighbors = 20, n_pcs=40)
sc.tl.umap(muscle)


# In[85]:


sc.pl.umap(muscle, color = ["factor.8", "factor.7", "Hours"], add_outline = True, size = 200)


# In[86]:


muscle.obs["Hours"] = muscle.obs["Hours"].astype("category")


# In[88]:


adata = muscle

import numpy as np
import matplotlib.pyplot as plt

fig, (ax1, ax2, ax3, ax4) = plt.subplots(1,4, figsize = (9,3))

info = adata.uns["dict_res_prop"]["factor.1"]
langs = np.array(info[2]).astype("str")
students = (np.array(info[0]).astype("int"))/(np.array(info[1]).astype("int")) * 100
ax1.bar(langs,students, color = "black")
ax1.set_xlabel("Hours", fontsize = "x-large")
ax1.set_ylabel("% activated cells", fontsize = "xx-large", labelpad=10)
ax1.set_title("Factor loading 1", fontsize = "x-large", pad=10)
ax1.set_ylim((0,35))
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(True)
#ax.spines['left'].set_visible(False)


info = adata.uns["dict_res_prop"]["factor.5"]
langs = np.array(info[2]).astype("str")
students = (np.array(info[0]).astype("int"))/(np.array(info[1]).astype("int")) * 100
ax2.bar(langs,students, color = "black")
ax2.set_xlabel("Hours", fontsize = "x-large")
#ax2.set_ylabel("% cells with gene set activation", fontsize = "x-large")
ax2.set_title("Factor loading 5", fontsize = "x-large", pad=10)
ax2.set_ylim((0,35))
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(True)

info = adata.uns["dict_res_prop"]["factor.54"]
langs = np.array(info[2]).astype("str")
students = (np.array(info[0]).astype("int"))/(np.array(info[1]).astype("int")) * 100
ax3.bar(langs,students, color = "black")
ax3.set_xlabel("Hours", fontsize = "x-large")
#ax3.set_ylabel("% cells with gene set activation", fontsize = "x-large")
ax3.set_title("Factor loading 54", fontsize = "x-large", pad=10)
ax3.set_ylim((0,35))
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(True)

info = adata.uns["dict_res_prop"]["factor.0"]
langs = np.array(info[2]).astype("str")
students = (np.array(info[0]).astype("int"))/(np.array(info[1]).astype("int")) * 100
ax4.bar(langs,students, color = "black")
ax4.set_xlabel("Hours", fontsize = "x-large")
#ax3.set_ylabel("% cells with gene set activation", fontsize = "x-large")
ax4.set_title("Factor loading 0", fontsize = "x-large", pad=10)
ax4.set_ylim((0,35))
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['bottom'].set_visible(True)

fig.tight_layout(pad = 4.0)

plt.savefig("results/muscle_factors.png", dpi = 600)


# In[89]:


fig, (ax1, ax2, ax3, ax4) = plt.subplots(1,4, figsize = (9,6))

degenes = []
genesofinterest = get_df_for_factor(muscle, 1).sort_values("z_score", ascending = False).head(20).index.to_list()
loadingscores = get_df_for_factor(muscle, 1).sort_values("z_score", ascending = False).head(20).z_score.to_list()
color_list = []
genes_plotting = pd.DataFrame(zip(genesofinterest, loadingscores), index = genesofinterest).reindex(genesofinterest)

Product = genes_plotting[0]

for i in Product:
    if i in degenes:
        color_list.append("black")
    else:
        color_list.append("white")

Quantity = genes_plotting[1]

ax1.barh(Product, Quantity, color = color_list, edgecolor = "black")
ax1.invert_yaxis()

ax1.set_ylabel("Top genes in factor loading", fontsize = "xx-large", labelpad=10)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(False)

import matplotlib.ticker as ticker
ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax1.grid(axis='x')

#plt.ylabel('Gene')
ax1.set_xlabel('Z-scores')


degenes = []
genesofinterest = get_df_for_factor(muscle, 5).sort_values("z_score", ascending = False).head(20).index.to_list()
loadingscores = get_df_for_factor(muscle, 5).sort_values("z_score", ascending = False).head(20).z_score.to_list()
color_list = []
genes_plotting = pd.DataFrame(zip(genesofinterest, loadingscores), index = genesofinterest).reindex(genesofinterest)

Product = genes_plotting[0]

for i in Product:
    if i in degenes:
        color_list.append("black")
    else:
        color_list.append("white")

Quantity = genes_plotting[1]

ax2.barh(Product, Quantity, color = color_list, edgecolor = "black")
ax2.invert_yaxis()

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(True)
ax2.spines['left'].set_visible(False)

import matplotlib.ticker as ticker
ax2.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax2.grid(axis='x')

#plt.ylabel('Gene')
ax2.set_xlabel('Z-scores')

degenes = []
genesofinterest = get_df_for_factor(muscle, 54).sort_values("z_score", ascending = False).head(20).index.to_list()
loadingscores = get_df_for_factor(muscle, 54).sort_values("z_score", ascending = False).head(20).z_score.to_list()
color_list = []
genes_plotting = pd.DataFrame(zip(genesofinterest, loadingscores), index = genesofinterest).reindex(genesofinterest)

Product = genes_plotting[0]

for i in Product:
    if i in degenes:
        color_list.append("black")
    else:
        color_list.append("white")

Quantity = genes_plotting[1]

ax3.barh(Product, Quantity, color = color_list, edgecolor = "black")
ax3.invert_yaxis()

ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(True)
ax3.spines['left'].set_visible(False)

ax3.xaxis.set_major_locator(ticker.MultipleLocator(4))
ax3.grid(axis='x')

#plt.ylabel('Gene')
ax3.set_xlabel('Z-scores')

degenes = []
genesofinterest = get_df_for_factor(muscle, 0).sort_values("z_score", ascending = False).head(20).index.to_list()
loadingscores = get_df_for_factor(muscle, 0).sort_values("z_score", ascending = False).head(20).z_score.to_list()
color_list = []
genes_plotting = pd.DataFrame(zip(genesofinterest, loadingscores), index = genesofinterest).reindex(genesofinterest)

Product = genes_plotting[0]

for i in Product:
    if i in degenes:
        color_list.append("black")
    else:
        color_list.append("white")

Quantity = genes_plotting[1]

ax4.barh(Product, Quantity, color = color_list, edgecolor = "black")
ax4.invert_yaxis()

ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['bottom'].set_visible(True)
ax4.spines['left'].set_visible(False)

import matplotlib.ticker as ticker
ax4.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax4.grid(axis='x')

#plt.ylabel('Gene')
ax4.set_xlabel('Z-scores')

fig.tight_layout(pad = 4.0)
plt.savefig("results/muscle_genes.png", dpi = 300)


# In[ ]:




