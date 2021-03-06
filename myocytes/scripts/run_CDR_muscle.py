
import anndata as ad
from pycdr.pycdr import run_CDR_analysis
from pycdr.perm import calculate_enrichment
from enrichment_utils.ontology_analysis import analyse_adata

INPUT = snakemake.input[0]
INPUT_FILE_GENE2GO = snakemake.input[1]
INPUT_FILE_ONTOLOGY = snakemake.input[2]

OUTPUT = snakemake.output[0]

x = ad.read(INPUT)

run_CDR_analysis(x, "Hours")

analyse_adata(x, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human", threshold_pvalue = 0.05, ontology_subset = "BP", prop = False)
factor_list = [i for i in x.uns["factor_loadings"].keys()]

calculate_enrichment(x, "Hours", factor_list, 100, "gene_short_name", 0.05)

x.write(OUTPUT)
