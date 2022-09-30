#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import anndata as ad
import scanpy as sc
from pycdr.perm import calculate_enrichment

muscle = ad.read(snakemake.input[0])

calculate_enrichment(muscle, "Hours", ["factor.1", "factor.0", "factor.8", "factor.5", "factor.7", "factor.54", "factor.8"], 100, "features", 0.1, 42)

#sc.pp.scale(muscle, max_value=10)
sc.tl.pca(muscle, svd_solver='arpack')
sc.pp.neighbors(muscle, n_neighbors = 20, n_pcs=50)
sc.tl.umap(muscle, random_state = 7)


muscle.obs["Hours"] = muscle.obs["Hours"].astype("category")

sc.set_figure_params(scanpy=True, dpi=80, dpi_save=600, format='png')
sc.pl.umap(muscle, color = ["factor.8", "factor.7", "Hours"], add_outline = True, size = 200, save="muscle")

muscle.write(snakemake.output[0])
