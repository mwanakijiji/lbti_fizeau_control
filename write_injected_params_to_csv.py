#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This reads in FITS files and extracts injected parameters (OPD, tip, tilt, x, y)
# and writes them to a csv

# Created 2020 Jan 1 by E.S.


# In[1]:


import glob
from astropy.io import fits
import pandas as pd


# In[2]:


# trial1 synthetic dataset

file_list = glob.glob("*trial5*.fits")

# initialize, considering elements strings
df = pd.DataFrame('', index=range(0,len(file_list)), 
                  columns=['filename','opd','tip','tilt','x_shift','y_shift'])

for f in range(0,len(file_list)):

    image, hdr = fits.getdata(file_list[f],0,header=True)
    
    df.iloc[f]["filename"] = str(file_list[f])
    df.iloc[f]["opd"] = hdr["OPD_UM"]
    df.iloc[f]["tip"] = hdr["TIPY_MAS"]
    df.iloc[f]["tilt"] = hdr["TILTXMAS"]
    
    # some datasets have these too
    if ("X_SHIFT" in list(hdr.keys())):
        df.iloc[f]["x_shift"] = hdr["X_SHIFT"]
        df.iloc[f]["y_shift"] = hdr["Y_SHIFT"]
        
    
# write out
df.to_csv("junk2.csv", index=False)

