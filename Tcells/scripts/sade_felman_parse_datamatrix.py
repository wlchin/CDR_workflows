#!/usr/bin/env python
# coding: utf-8

# In[4]:


import anndata as ad
import pandas as pd
import glob
import os

humnum = pd.read_csv(snakemake.input[0])["Cell Name"]

# required to get the rows
unos = pd.read_csv(snakemake.input[2], 
            nrows = 0, skiprows=0,
            engine = "python", 
           sep = "\t",
           index_col=False)
listname = unos.columns.to_list()
listname.pop(0)


chunksize = 1000
with pd.read_csv(snakemake.input[1], chunksize=chunksize, sep = '\t', header = None, index_col=0) as reader:
    for i, chunk in enumerate(reader):
        a = chunk.drop(16292, axis = 1)
        a.index.name = "genes"
        a.astype('float64')
        a.columns = listname
        selecto = a.columns[a.columns.isin(humnum)]
        uiu = a[selecto]
        fn = "test" + str(i) + ".csv"
        uiu.to_csv(fn)


hum = glob.glob("*.csv")
new = pd.concat([pd.read_csv(i) for i in hum])
new = new.set_index("genes")

[os.remove(i) for i in hum]

new.to_pickle(snakemake.output[0])

