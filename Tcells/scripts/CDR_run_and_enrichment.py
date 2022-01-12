#!/usr/bin/env python
# coding: utf-8

import anndata as ad
from pycdr.pycdr import run_CDR_analysis
from enrichment_utils.ontology_analysis import analyse_adata

x = ad.read(snakemake.input[0])

run_CDR_analysis(x, "response_pheno")

analyse_adata(x, snakemake.input[1], 
                  snakemake.input[2], 
                  "human", 
                  threshold_pvalue = 0.05, 
                  ontology_subset = "BP", 
                  prop = False)

x.write(snakemake.output[0])

