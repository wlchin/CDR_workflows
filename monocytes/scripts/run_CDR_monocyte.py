
import anndata as ad
from pycdr.pycdr import run_CDR_analysis

INPUT = snakemake.input[0]

OUTPUT = snakemake.output[0]

mono = ad.read(INPUT)
run_CDR_analysis(mono, "stim")

mono.write(OUTPUT)
