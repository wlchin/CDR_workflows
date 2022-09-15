

import anndata as ad
from pycdr.utils import filter_genecounts_numcells, filter_genecounts_percent
import scanpy as sc

mono = ad.read(snakemake.input[0])

mono = filter_genecounts_percent(mono, 0.01, 1)
mono = filter_genecounts_numcells(mono, 0, 100)

mono.raw = mono.raw.to_adata()

sc.pp.log1p(mono)
sc.pp.scale(mono, zero_center = False)

mono.write(snakemake.output[0])

