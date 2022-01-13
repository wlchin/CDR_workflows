#!/usr/bin/env python
# coding: utf-8

import anndata as ad
import scanpy as sc
from pycdr.utils import filter_genecounts_percent
import numpy as np

tirosh = ad.read(snakemake.input[0])
sade_felman = ad.read(snakemake.input[1])

sc.pp.log1p(sade_felman)
sc.pp.log1p(tirosh)

sc.pp.highly_variable_genes(sade_felman)
sc.pp.highly_variable_genes(tirosh)

sade_felman = sade_felman[:, sade_felman.var.highly_variable]
tirosh = tirosh[:, tirosh.var.highly_variable]

sc.pp.scale(sade_felman, zero_center = False)
sc.pp.scale(tirosh, zero_center = False)

combined_dataset = sade_felman.concatenate(tirosh)

combined_dataset.obs.dropna(axis = 1, inplace=True)

combined_dataset.var["features"] = combined_dataset.var["features-1"]

combined_dataset_filtered = filter_genecounts_percent(combined_dataset, 0.01, 1)
combined_dataset_filtered = combined_dataset_filtered.copy()

print("dimensions: ", str(combined_dataset_filtered.shape))
combined_dataset_filtered.write(snakemake.output[0])




