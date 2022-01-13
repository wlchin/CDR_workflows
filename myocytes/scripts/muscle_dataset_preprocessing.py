import scanpy as sc
import anndata as ad
from pycdr.utils import filter_genecounts_numcells, filter_genecounts_percent

muscle = ad.read(snakemake.input[0])
muscle.var_names_make_unique()
muscle.var["features"] = muscle.var["gene_short_name"]
filtered_by_numcells = filter_genecounts_numcells(muscle, 1, 50)

sc.pp.log1p(filtered_by_numcells)
sc.pp.scale(filtered_by_numcells, zero_center = False)

filtered_by_numcells.write(snakemake.output[0])




