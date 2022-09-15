import anndata as ad
import scanpy as sc
import pandas as pd
from enrichment_utils.ontology_analysis import analyse_list

INPUT_FILE = snakemake.input[0]
INPUT_FILE_GENE2GO = snakemake.input[1]
INPUT_FILE_ONTOLOGY = snakemake.input[2]

OUTPUT_DE_ENRICHMENT_LOW = snakemake.output[0]
OUTPUT_DE_ENRICHMENT_HIGH = snakemake.output[1]
OUTPUT_DE_LIST_LOW = snakemake.output[2]
OUTPUT_DE_LIST_HIGH = snakemake.output[3]


def get_DE_genes(adata, factor, pheno1, pheno2, lfc, fdr):
    adata.var_names_make_unique()
    adata_filtered = adata[(adata.obs[factor] == pheno1) | (adata.obs[factor] == pheno2)]
    adata_filtered.obs[factor] = adata_filtered.obs[factor].astype("string")

    sc.tl.rank_genes_groups(adata_filtered, factor, method='wilcoxon', key_added = "wilcoxon")

    wcup = sc.get.rank_genes_groups_df(adata_filtered, 
                                       group=pheno1, key='wilcoxon', 
                                       pval_cutoff=fdr, log2fc_min=lfc)['names']

    wcdown = sc.get.rank_genes_groups_df(adata_filtered, 
                                         group=pheno1, key='wilcoxon', 
                                         pval_cutoff=fdr, log2fc_max = -lfc)['names']

    upgenes = wcup.to_list()
    downgenes = wcdown.to_list()
    comb = upgenes + downgenes

    return comb, upgenes, downgenes

def combine_output_of_lists(set1, set2, outputfile):
    import pandas as pd
    set1 = set(set1)
    set2 = set(set2)
    listo = set1.union(set2)
    x = pd.DataFrame(listo)
    x.columns = ["Terms"]
    x.to_csv(outputfile)


x = ad.read(INPUT_FILE)
adata = x.raw.to_adata() # remember that 
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
x.raw = adata

a,b,c = get_DE_genes(x, "stim", "CTRL", "STIM", 0.5, 0.1) # low log
e,f,g = get_DE_genes(x, "stim", "CTRL", "STIM", 1.0, 0.05) # strict

gos, enrichment1_low = analyse_list(a, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human")
gos, enrichment2_low = analyse_list(b, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human")
gos, enrichment3_high = analyse_list(e, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human")
gos, enrichment4_high = analyse_list(f, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human")

combine_output_of_lists(enrichment1_low, enrichment2_low, OUTPUT_DE_ENRICHMENT_LOW)
combine_output_of_lists(enrichment3_high, enrichment4_high, OUTPUT_DE_ENRICHMENT_HIGH)

high_thres = pd.DataFrame(g)
high_thres.to_csv(OUTPUT_DE_LIST_HIGH)

low_thres = pd.DataFrame(c)
low_thres.to_csv(OUTPUT_DE_LIST_LOW)




