
# Spectral detection of condition-specific biological pathways in single-cell gene expression data.

This repository contains the snakemake workflows required to generate the results in the CDR-g manuscript. 

These workflows contain steps that are containerised. Therefore, involve the --run-singularity flag when runnin snakemake. 

# Recommended dependencies:

Scanpy, bbknn and enrichment_utils are required for data proprocessing and visualisation of results. Please refer to the [CDR-g page](https://github.com/wlchin/pycdr) for additional installation information.

# Functional annotation - required file versions

Enrichment analysis in these workflows depends on the PANTHER GO-ontology (version 16.0) and the NCBI gene2go (dated 2022-01-03). These files are provided in this repository.  