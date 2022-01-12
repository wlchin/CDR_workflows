

import anndata as ad
from pycdr.utils import filter_genecounts_numcells, filter_genecounts_percent
import scanpy as sc

mononess = ad.read(snakemake.input[0])

duo = filter_genecounts_percent(mononess, 0.01, 1)
uno = filter_genecounts_numcells(duo, 0, 100)

uno.raw = uno.raw.to_adata()

sc.pp.log1p(uno)
sc.pp.scale(uno, zero_center = False)

uno.write(snakemake.output[0])

