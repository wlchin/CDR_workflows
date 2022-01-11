
import anndata as ad
from pycdr.pycdr import run_CDR_analysis
from pycdr.perm import calculate_enrichment
from enrich_cdr.encdr import do_enrichment_all

INPUT = snakemake.input[0]
INPUT_FILE_GENE2GO = snakemake.input[1]
INPUT_FILE_ONTOLOGY = snakemake.input[2]

OUTPUT = snakemake.output[0]

mono = ad.read(INPUT)
run_CDR_analysis(mono, "stim")
do_enrichment_all(mono, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human", threshold_pvalue = 0.05, ontology_subset = "BP", prop = False)
factor_list = [i for i in mono.uns["factor_loadings"].keys()]
calculate_enrichment(mono, "stim", factor_list, 100, "features", 0.05)

mono.write(OUTPUT)
