#!/usr/bin/env python
# coding: utf-8

import anndata as ad
import pandas as pd
import glob
import os

# get list of CD8 T cells
CD8Tcells = pd.read_csv(snakemake.input[0])["Cell Name"]

# parse raw file to get colnames
rawfile = pd.read_csv(snakemake.input[2], 
            nrows = 0, skiprows=0,
            engine = "python", 
           sep = "\t",
           index_col=False)
colnames_of_raw_datafile = rawfile.columns.to_list()
colnames_of_raw_datafile.pop(0)


chunksize = 1000
with pd.read_csv(snakemake.input[1], chunksize=chunksize, sep = '\t', header = None, index_col=0) as reader:
    for i, chunk in enumerate(reader):
        a = chunk.drop(16292, axis = 1)
        a.index.name = "genes"
        a.astype('float64')
        a.columns = colnames_of_raw_datafile
        select_relevant_cells = a.columns[a.columns.isin(CD8Tcells)]
        processed_chunked = a[select_relevant_cells]
        chunk_filenamefilename = "chunk" + str(i) + ".csv"
        processed_chunked.to_csv(chunk_filename)

csv_chunks = glob.glob("*.csv")
new = pd.concat([pd.read_csv(i) for i in csv_chunks])
new = new.set_index("genes")


[os.remove(i) for i in csv_chunks]
new.to_pickle(snakemake.output[0])

