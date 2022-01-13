#!/usr/bin/env python
# coding: utf-8

# # No prior information experiments (Fig 2 + S5)
# 
# We use Benchmark 2 (derived from recently deposited RNA structures in the PDB) to assess prediction performance when prior experimental information is ignored.
# 
# Pdbs are available at https://purl.stanford.edu/sq987cc0358, with one folder per RNA.

# In[1]:


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# In[2]:


# No bootstrapping.
benchmark2_nobootstrap = pd.read_csv('benchmark2_nobootstrap.csv')
# Bootstrapping.
benchmark2_bootstrap = pd.read_csv('benchmark2_bootstrap.csv')


# # RMSD of best-scoring structural model, bootstrapped (Fig 2D)

# In[3]:


plt.figure(figsize=(12, 4))
fig = plt.gcf()
fig.patch.set_alpha(0)

plt.errorbar(x=benchmark2_bootstrap.year, 
             y=benchmark2_bootstrap.bo1_med, 
             yerr=benchmark2_bootstrap[['bo1_low_err', 'bo1_high_err']].T.values, 
             ls='', marker='o', markersize=4, capsize=3, color='r', ecolor='r')
plt.ylim([10, 19])
#plt.title('Best overall performance, new benchmark')

plt.ylabel('RMSD')
plt.xlabel('Year scoring function published')
xmin = 2006
xmax = 2022
plt.xlim([xmin, xmax])
plt.grid(True, lw=0.25)
plt.tight_layout()
plt.show()


# # RMSD of best-scoring structural model, no bootstrapping (Fig S5A)

# In[4]:


plt.figure(figsize=(12, 4))
fig = plt.gcf()
fig.patch.set_alpha(0)

plt.errorbar(x=benchmark2_nobootstrap.year, 
             y=benchmark2_nobootstrap.bo1_med, 
             ls='', marker='+', markersize=10, capsize=3, color='r', ecolor='r')
plt.ylim([10, 19])
#plt.title('Best overall performance, new benchmark')

plt.ylabel('RMSD')
plt.xlabel('Year scoring function published')
xmin = 2006
xmax = 2022
plt.xlim([xmin, xmax])
plt.grid(True, lw=0.25)
plt.tight_layout()
plt.show()


# # Lowest RMSD among 10 best-scoring structural models, bootstrapped (Fig S5B)

# In[5]:


plt.figure(figsize=(12, 4))
fig = plt.gcf()
fig.patch.set_alpha(0)

plt.errorbar(x=benchmark2_bootstrap.year, 
             y=benchmark2_bootstrap.bo10_med, 
             yerr=benchmark2_bootstrap[['bo10_low_err', 'bo10_high_err']].T.values, 
             ls='', marker='o', markersize=4, capsize=3, color='r', ecolor='r')
plt.ylim([8, 14.5])
#plt.title('Best performance in top 10, new benchmark')

plt.ylabel('RMSD')
plt.xlabel('Year scoring function published')
xmin = 2006
xmax = 2022
plt.xlim([xmin, xmax])
plt.grid(True, lw=0.25)
plt.tight_layout()
plt.show()


# # Lowest RMSD among 10 best-scoring structural models, no bootstrapping (Fig S5C)

# In[6]:


plt.figure(figsize=(12, 4))
fig = plt.gcf()
fig.patch.set_alpha(0)

plt.errorbar(x=benchmark2_nobootstrap.year, 
             y=benchmark2_nobootstrap.bo10_med, 
             ls='', marker='+', markersize=10, capsize=3, color='r', ecolor='r')
plt.ylim([8, 14.5])
#plt.title('Best performance in top 10, new benchmark')

plt.ylabel('RMSD')
plt.xlabel('Year scoring function published')
xmin = 2006
xmax = 2022
plt.xlim([xmin, xmax])
plt.grid(True, lw=0.25)
plt.tight_layout()
plt.show()


# In[ ]:




