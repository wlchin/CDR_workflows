
import anndata as ad
from pycdr.pycdr import run_CDR_analysis
from pycdr.perm import calculate_enrichment
from enrichment_utils.ontology_analysis import analyse_adata

INPUT = snakemake.input[0]
INPUT_FILE_GENE2GO = snakemake.input[1]
INPUT_FILE_ONTOLOGY = snakemake.input[2]

OUTPUT_enrichment = snakemake.output[0]
OUTPUT_enrichment_proportions = snakemake.output[1]

mono = ad.read(INPUT)
#run_CDR_analysis(mono, "stim")
analyse_adata(mono, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human", threshold_pvalue = 0.05, ontology_subset = "BP", prop = False)

mono.write(OUTPUT_enrichment)

factor_list = [i for i in mono.uns["factor_loadings"].keys()]
calculate_enrichment(mono, "stim", factor_list, 100, "features", 0.05)

mono.write(OUTPUT_enrichment_proportions)
