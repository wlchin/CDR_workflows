#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import pandas as pd

phenodata = pd.read_csv(snakemake.input[0], sep = "\t", header = 0, index_col=0)
clean = phenodata[phenodata.columns[0:6]]
clean.index = clean.title

otheraddtional = clean["characteristics: patinet ID (Pre=baseline; Post= on treatment)"].str.split("_", expand = True)
clean["stage"] = otheraddtional[0]
clean["res"] = clean["characteristics: response"] + "_" + clean["stage"] 

# create anndata from pickled df
genedf = pd.read_pickle(snakemake.input[1])
sade_felman = sc.AnnData(genedf)
sade_felman = sade_felman.T

sade_felman.obs["title"] = sade_felman.obs.index
pheno_reduced = clean[clean.title.isin(sade_felman.obs.index)]
pheno_reduced.index.name = "cell"
sade_felman.obs = sade_felman.obs.merge(pheno_reduced)

sade_felman.obs.index = sade_felman.obs.title
sade_felman.obs.index.name = "cellnames"

df = clean[clean.index.isin(sade_felman.obs.index)]
sade_felman.obs["response_pheno"] = df["res"]

sade_felman.write(snakemake.output[0])
sade_felman.obs.to_csv(snakemake.output[1])

