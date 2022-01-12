import anndata as ad
import pandas as pd

INPUT_CDR = snakemake.input[0]
INPUT_DE_LIST = snakemake.input[1]

OUTPUT_CDR_ONLY = snakemake.output[1]
OUTPUT_DE_ONLY = snakemake.output[2]
OUTPUT_INTERSECTION_ONLY = snakemake.output[0]


def retrieve_terms_for_selected_factors(x, factorlist):
    # earmarked removal move to other package
    en = x.uns["enrichment_results"]

    terms = []
    

    for i in factorlist:
        termie = en[i]
        terms.append(termie)

    flat_list = set([item for sublist in terms for item in sublist])

    len(set(flat_list))
    
    return set(flat_list) 

CDR_obj = ad.read(INPUT_CDR)

factor_list = [i for i in CDR_obj.uns["factor_loadings"].keys()]
CDR_set = retrieve_terms_for_selected_factors(CDR_obj, factor_list)

other = pd.read_csv(INPUT_DE_LIST)
DEset = set(other.Terms)

len(CDR_set), len(DEset - CDR_set), len(CDR_set - DEset)

pd.DataFrame(CDR_set - DEset).to_csv(OUTPUT_CDR_ONLY)
pd.DataFrame(DEset - CDR_set).to_csv(OUTPUT_DE_ONLY)
pd.DataFrame(DEset.intersection(CDR_set)).to_csv(OUTPUT_INTERSECTION_ONLY)




