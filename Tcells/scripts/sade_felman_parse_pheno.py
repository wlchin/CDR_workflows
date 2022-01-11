#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc
import pandas as pd


# In[1]:


phenodata = pd.read_csv(snakemake.input[0], sep = "\t", header = 0, index_col=0)
cleaner = phenodata[phenodata.columns[0:6]]
cleaner.index = cleaner.title
otheraddtional = cleaner["characteristics: patinet ID (Pre=baseline; Post= on treatment)"].str.split("_", expand = True)
cleaner["stage"] = otheraddtional[0]
cleaner["totality_as"] = cleaner["characteristics: response"] + "_" + cleaner["stage"]


# In[4]:


cleaner.head()


# In[5]:


cleaner["totality_as"].value_counts()


# In[ ]:


## get the gene DF

genedf = pd.read_pickle(snakemake.input[1])


# In[8]:


sade_felman = sc.AnnData(genedf)
sade_felman = sade_felman.T


sade_felman.obs["title"] = sade_felman.obs.index
smolpheno = cleaner[cleaner.title.isin(sade_felman.obs.index)]
smolpheno.index.name = "cell"
sade_felman.obs = sade_felman.obs.merge(smolpheno)


# In[9]:


sade_felman.obs.index = sade_felman.obs.title
sade_felman.obs.index.name = "cellnames"


# In[10]:


totalityas = cleaner[cleaner.index.isin(sade_felman.obs.index)]
sade_felman.obs["response_pheno"] = totalityas["totality_as"]
sade_felman.obs.response_pheno.value_counts()


# In[ ]:


sade_felman.write(snakemake.output[0])
sade_felman.obs.to_csv(snakemake.output[1])

