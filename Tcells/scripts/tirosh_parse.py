#!/usr/bin/env python
# coding: utf-8

# In[2]:


import anndata as ad
import pandas as pd


# In[3]:


tpm = ad.read_csv(snakemake.input[0])


# In[4]:


onlyyou = tpm.T


# In[5]:


#onlyyou.write("cleantirosh.h5ad")


# In[6]:


hum = pd.read_csv(snakemake.input[1])


# In[7]:


hum.index = hum.cells


# In[8]:


onlyyou.obs = hum


# In[9]:


onlyyou.obs.index.name = "features"


# In[10]:


onlyyou.var["features"] = onlyyou.var.index


# In[11]:


onlyyou.obs["cell.types"].unique()


# In[12]:


CD8only = onlyyou[onlyyou.obs["cell.types"] == "T.CD8"]


# In[14]:


CD8only.obs["response_pheno"] = CD8only.obs["treatment.group"]


# In[16]:


CD8only.write(snakemake.output[0])


# In[17]:


CD8only

