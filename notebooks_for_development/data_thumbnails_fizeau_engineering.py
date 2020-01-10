#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This notebook makes thumbnails of lots of Fizeau frames so as to isolate what
# I want to do engineering tests

# created from a parent, 2020 Jan. 9 by E.S.


# In[4]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os.path
import glob


# In[2]:


from lmircam_tools import *
from lmircam_tools import overlap_psfs


# In[12]:


# get all filenames in the directory

dir_name = ('.')
file_names_all = glob.glob(dir_name + "*.fits")


# In[13]:


file_names_all


# In[22]:


for f in range(0,len(file_names_all)): # full dataset in HOSTS footprint: (4249,11497)
    
    image, header = fits.getdata(file_names_all[f],0,header=True) 
    
    # a quick and dirty background subtraction and bad pixel fixing
    #TBD
    
    # make plot
    plt.imshow(image, origin="lower", cmap="gray") # PSF cut-out
    plt.ylabel('y')
    plt.xlabel('x')
    plt.suptitle(str(os.path.basename(file_names_all[f])))

    plt.tight_layout()
    plt.savefig("thumbnail_"+str(os.path.basename(file_names_all[f]).split(".")[0])+".jpeg", dpi=100, overwrite=False)
    plt.clf()
    
    print('Frame '+str(os.path.basename(file_names_all[f]))+' done...')
    

