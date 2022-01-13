#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import anndata as ad
import scanpy as sc
from pycdr.perm import calculate_enrichment

muscle = ad.read(snakemake.input[0])

calculate_enrichment(muscle, "Hours", ["factor.1", "factor.0", "factor.8", "factor.5", "factor.7", "factor.54", "factor.8"], 100, "features", 0.1, 19)

sc.pp.scale(muscle, max_value=10)
sc.tl.pca(muscle, svd_solver='arpack')
sc.pp.neighbors(muscle, n_neighbors = 20, n_pcs=40)
sc.tl.umap(muscle)

sc.pl.umap(muscle, color = ["factor.8", "factor.7", "Hours"], add_outline = True, size = 200, save = ".png")

muscle.obs["Hours"] = muscle.obs["Hours"].astype("category")

muscle.write(snakemake.output[0])
