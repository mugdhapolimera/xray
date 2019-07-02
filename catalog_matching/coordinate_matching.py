#!/usr/bin/env python
# coding: utf-8

# # Cross-Matching VOICE-VST Photometric and Spectroscopic Redshift Catalogs
# 
# Zack Hutchens - June 2019

# In[1]:


import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.stats import mode, ttest_rel, linregress
from astropy.io import fits
import os

#os.chdir(r'C:\Users\mugdhapolimera\github\xray\catalog_matching')
def allUnique(x):
    seen = set()
    return not any(i in seen or seen.add(i) for i in x)


# In[54]:


spec_file = r"../../SDSS_spectra/ECO+RESOLVE_filter_new.csv"
#phot_file = r"~/Desktop/3XMM_DR8_slim.csv"

# Read in RA + Dec from both catalogs
sdf = pd.read_csv(spec_file)
sdf
#pdf = np.genfromtxt(phot_file, delimiter = ' ')


# In[111]:


phot_file = r'../../SDSS_spectra/Broadlineagn.fits'
#pdf = pd.read_csv(phot_file, delimiter = ',')
from astropy.io import fits
from astropy.table import Table
pdf_fits = Table.read(phot_file)
pdf = pdf_fits.to_pandas()
pdf = pdf.rename(columns = {'DESIGNATION':'Name',"RA":'radeg', 'DEC':'dedeg' })
pdf.Name = [x.decode('utf-8') for x in pdf.Name]
pdf


# In[100]:


#indices = np.where(np.array(pdf["bestEstimate"]) >= 0)

# in the pdf frame, there are some z = -99. Remove these using [indices].
sra = np.array(sdf["radeg"])
sdec = np.array(sdf["dedeg"])
pra = np.array(pdf["radeg"])#[indices]
pdec = np.array(pdf["dedeg"])#[indices]
len(sra), len(pra)


# In[101]:


scatalog = SkyCoord(ra=sra*u.degree, dec=sdec*u.degree)
pcatalog = SkyCoord(ra=pra*u.degree, dec=pdec*u.degree)
idx, d2d, d3d = match_coordinates_sky(scatalog, pcatalog, nthneighbor=1)

matches = pcatalog[idx]


# Now, `idx` contains the indices in `pcatalog` which are closest matches to the targets in `scatalog.`

# In[102]:


matches, len(matches)


# In[103]:


len(idx), max(idx), min(idx)


# ## Cut out only those that are the closest matches: `idx` is not unique

# In[104]:


def find_closest_matches(scatalog, pcatalog, idx, d2d):
    """
    function `find_closest_matches` --- Z. Hutchens, June 2019.
    Find the exact closest match for every target in the matched results of astropy.match_coordinates_*.
    Arguments: 
        scatalog: smaller catalog (astropy.SkyCoord instance)
        pcatalog: larger catalog (astropy.SkyCoord instance) to which `scatalog` has been matched.
        idx: the indices in `pcatalog` that contains the closest matches for each target in `scatalog`. (np.array)
        d2d: the on-sky distances between matched targets from astropy.match_coordinates_*. Shape matches idx. (np.array)
    Returns:
        data frame containing all targets in pcatalog and exact-closest matching targets in scatalog.
    """
    newdf = []
    for place in range(0, len(pcatalog.ra.value)):
        ind_in_scatalog = np.where(np.array(idx)==place)
            
        if len(ind_in_scatalog[0]) == 0:
            newdf.append([pcatalog.ra.value[place], pcatalog.dec.value[place], -99, -99, -99])
                          
                          #, pz[place],
                          #   -99, np.nan, np.nan, np.nan, np.nan,np.nan, np.nan, np.nan, np.nan, np.nan,np.nan,
                         #np.nan])
        else:
            dist = d2d.arcsec[ind_in_scatalog] # this has index values we don't care about, need orig indices
                
            # Combine dist, ind_in_scatalog into tuple-filled array
            dist_with_orig_ind = [(i,j) for i,j in zip(ind_in_scatalog[0], dist)]
                
            # find the minimum separation in the array
            minv = dist_with_orig_ind[0]
            for tup in dist_with_orig_ind:
                if tup[1] < minv[1]:
                    minv = tup
                
            index = minv[0]
            # now we have the (index, sep) of the smallest separation between matched coordinates at index
            # minv[0] in scatalog.
            # Now, append the items of interest to the new dataframe.
                
            newdf.append([sdf.name[index], scatalog.ra.value[index], scatalog.dec.value[index],
                          pdf.Name[place], pcatalog.ra.value[place], pcatalog.dec.value[place], d2d[index].arcsec]) 
                          
                          #pz[place], sz[index],
                         #sdss_mags["u"][index], sdss_mags["g"][index], sdss_mags["r"][index], 
                         # sdss_mags["i"][index], sdss_mags["z"][index], bes_mags["U"][index], bes_mags["B"][index],
                         # bes_mags["V"][index], bes_mags["R"][index], bes_mags["I"][index], s_names[index]])
    
    #cnames = ["photRA", "photDec", "specRA", "specDec", "separation", "photZ", "specZ", "sdss_absmagu",
    #          "sdss_absmagg", "sdss_absmagr", "sdss_absmagi", "sdss_absmagz", 
    #          "bessell_absmagU", "bessell_absmagB", "bessell_absmagV", "bessell_absmagR", "bessell_absmagI",
    #         "specName"]
    cnames = ["eco+res_name", "radeg", "dedeg","agn_name", "agn_ra", "agn_dec", "separation"]
    return pd.DataFrame(newdf, columns = cnames)


# In[112]:


matched = find_closest_matches(scatalog, pcatalog, idx, d2d)


# In[113]:


print (len(matched[matched.agn_ra != -99.0][matched.separation < 7]))
matched = matched[matched.agn_ra != -99.0][matched.separation < 7]
matched


# In[114]:


matched.to_csv("/afs/cas.unc.edu/users/m/u/mugpol/github/xray/catalog_matching/BroadlineAgnMatched.csv", index=False)

