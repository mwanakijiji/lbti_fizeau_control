#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This scrapes data from the 2018 Fizeau AC Her data, 2018 July 9 by E.S.


# In[1]:


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
import asciitable
import sys
import glob, os


# In[21]:


# read in a file and extract a meta-data field from the FITS header

dirTreeStem = ('./')
fileNameStem = ('lm_190224_')

os.chdir(dirTreeStem)

counter = 0

for f in glob.glob("*.fits"): # loop over filenames
    counter += 1 # count loops
    
    counterStart = 60
    counterStop = 11500
    if (np.logical_or(counter<counterStart,counter>counterStop)):
        continue
    
    if (counter%50 == 0):
        print('Reading in header info from frame '+str(f)+'...')
    
    # nice continuous printing of status, but can get stuck 
    #print('\r', 'Reading in header info from frame '+str(f)+'...', end='')
    #sys.stdout.flush()
    
    image, header = fits.getdata(dirTreeStem+f,0,header=True)
    
    
    # initialize dictionary with empty lists for each key, if this is first loop iteration
    if (counter==counterStart):
        data_ascii = {}
        for k in header:
            data_ascii[k] = [] 
        data_ascii['NODPOS'] = [] # manually add in header keys that don't appear in all frames
        data_ascii['DATAFLAG'] = []
        data_ascii['BADROWS'] = []
        data_ascii['FRAMENUM'] = []
    

    # append data from frames   
    for l in data_ascii:
        if l in header:
            data_ascii[l].append(header[l])
        else:
            if (l!='FRAMENUM'):
                data_ascii[l].append(np.nan) # if a flag doesn't appear in the header
            else:
                split_underscore = f.split("_") # to get frame number
                #print(split_underscore)
                split_final = split_underscore[2].split(".")
                data_ascii[l].append(split_final[0])
                
    #print(np.shape(data_ascii))


# In[22]:


# convert to pandas dataframe

data_ascii_df = pd.DataFrame.from_dict(data_ascii, orient='columns')


# In[24]:


# write to file

data_ascii_df.to_csv('junk.csv')

