#!/usr/bin/env python
# coding: utf-8

import anndata as ad
import scanpy as sc
import numpy as np
import bbknn
from pycdr.perm import calculate_enrichment


sc.set_figure_params(dpi_save = 300)

adata_all = ad.read(snakemake.input[0])

calculate_enrichment(adata_all, "response_pheno", 
                     ["factor.169", "factor.1192", "factor.10", "factor.13"], 
                     100, "features", 0.10, 19)


sc.tl.pca(adata_all)
bbknn.bbknn(adata_all, batch_key = "batch")
sc.tl.leiden(adata_all)
bbknn.ridge_regression(adata_all, batch_key=['batch'], confounder_key=['leiden'])

sc.tl.pca(adata_all)
sc.pp.neighbors(adata_all, n_neighbors=30, n_pcs=25)
sc.tl.leiden(adata_all, resolution = 1.0)
sc.tl.tsne(adata_all)

sc.pl.tsne(adata_all, color=["leiden", "factor.1192", "MT1E"], legend_loc='on data', save = "_MT1E.png")


cl1 = adata_all[adata_all.obs['leiden'].isin(["5"]),:]


sc.tl.rank_genes_groups(cl1, "factor.1192", method="wilcoxon")
sc.pl.rank_genes_groups(cl1, n_genes=25, sharey=False, save="_DEresult.png")


adata_all.obs["Exhaustion"] = adata_all.obs["factor.13"]
adata_all.obs["Cell division"] = adata_all.obs["factor.10"]
adata_all.obs["IFN response"] = adata_all.obs["factor.169"]
adata_all.obs["Hypoxia"] = adata_all.obs["factor.1192"]

sc.pl.tsne(adata_all, color = ["Exhaustion", "TIGIT", "PDCD1", 
                               "Cell division", "CDK1", "CDCA8", 
                               "IFN response", "MX1", "IFI44L", 
                               "Hypoxia", "BAD", "CCNB1"], 
           ncols=3, wspace = 0.3, save = "_genesets.png", legend_loc = None)





sc.tl.leiden(adata_all, resolution = 5.0) # 3 and 15 for this.
sc.pl.tsne(adata_all, color=["leiden"], save="leidenbignonum.png")

cl1 = adata_all[adata_all.obs['leiden'].isin(["1", "6", "19"]),:]

sc.tl.rank_genes_groups(cl1, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(cl1, n_genes=25, sharey=False, save="_DEresult_multicluster.png")