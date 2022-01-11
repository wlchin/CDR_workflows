#!/usr/bin/env python
# coding: utf-8

# In[98]:


import anndata as ad
import scanpy as sc
from pycdr.pp import filter_genecounts_percent, filter_genecounts_numcells
import numpy as np


# In[97]:


tirosh = ad.read(snakemake.input[0])
sade_felman = ad.read(snakemake.input[1])


# In[99]:


sc.pp.log1p(sade_felman)
sc.pp.log1p(tirosh)


# In[63]:


tirosh.obs.head() # note that this requires a response pheno.


# In[65]:


sc.pp.highly_variable_genes(sade_felman)
sc.pp.highly_variable_genes(tirosh)


# In[66]:


sade_felman = sade_felman[:, sade_felman.var.highly_variable]
tirosh = tirosh[:, tirosh.var.highly_variable]


# In[90]:


sc.pp.scale(sade_felman, zero_center = False)
sc.pp.scale(tirosh, zero_center = False)


# In[67]:


totality = sade_felman.concatenate(tirosh)


# In[68]:


totality.obs.dropna(axis = 1, inplace=True)


# In[85]:


totality


# In[86]:


totality.var["features"] = totality.var["features-1"]


# In[93]:


print("TCF7" in totality.var.index)
print("HAVCR2" in totality.var.index)
print("CTLA4" in totality.var.index)
print("TIGIT" in totality.var.index)


# In[184]:


13536 * 0.01


# In[89]:


a = filter_genecounts_percent(totality, 0.01, 1)


# In[95]:



print("TCF7" in a.var.index)
print("HAVCR2" in a.var.index)
print("CTLA4" in a.var.index)
print("TIGIT" in a.var.index)
print("PDCD1" in a.var.index)


# In[96]:

a = a.copy()

print("this is the result of filtering")
print(a.shape)

a.write(snakemake.output[0])


# In[ ]:




