#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Made for matching catalogs.
"""
import numpy as np
from astropy import units as u
import pandas as pd
import matplotlib.pyplot as plt
import os
arcsec = 1/3600.


def catmatch_act(ra1, dec1, ra2, dec2, rfluxes, fullflux, goodm1):
    '''
    ra2 to match to (GSW)
    ra1 to be matched (3XMM)
    '''
    match_7 = []
    matchedind = []
    matchra_diff = []
    matchdec_diff = []
    matchxflux =[]
    matchrflux = []
    match_7_dists = []
    for i in range(len(ra1)): #np.int64(goodm1): #
        racop = np.copy(ra2)
        #if i%1000 == 0:
        #    print(i)
        if ra1[i] >(360-7*arcsec):
            wrapra2 = np.where(racop<7*arcsec)
            racop[wrapra2] +=360
        elif ra1[i] <7*arcsec:
            wrapra2 = np.where(racop>(360-7*arcsec))
            racop[wrapra2] -=360

        delta_ra = (ra1[i] - racop)*np.cos(np.radians(dec1[i]))
        delta_dec = (dec1[i] - dec2)
        delta_rad = np.sqrt(delta_ra**2 + delta_dec**2)
        good_gals = np.where(delta_rad < 7*arcsec)[0]
        n_cands = len(good_gals)
        fullfl = fullflux[i]
        if n_cands >0:
            print(delta_ra[good_gals]/arcsec, delta_dec[good_gals]/arcsec, delta_rad[good_gals]/arcsec)
        if n_cands > 1:
            #pick the source with the greatest r-band flux
            comprflux = []
            for gal in good_gals:
                rflux = rfluxes[gal]
                comprflux.append(rflux)
            comprflux = np.array(comprflux)
#            print np.where(comprflux == np.max(comprflux))
            brightest = np.where(comprflux == np.min(comprflux))[0][0] #changed to min because using absmag
            galrflux = comprflux[brightest]
            good_gals = good_gals[brightest]
        elif len(good_gals) == 1:
            galrflux = rfluxes[good_gals]
            good_gals = good_gals[0]
        else:
            continue
        delta_ra_good = delta_ra[good_gals]
        delta_dec_good = delta_dec[good_gals]
        alreadymatched = np.where(match_7 == good_gals)[0]
        currentdist = delta_rad[good_gals]
        #print(delta_ra_good/arcsec, delta_dec_good/arcsec, currentdist/arcsec)
        if n_cands >1 and len(alreadymatched) >0:
            print('both messy')
        if len(alreadymatched) >0:
            #print(alreadymatched)
            #if an X-ray source has already been matched
            # take the one with largest full flux
            #
            alreadymatched= np.int64(alreadymatched[0])
            alreadymatchedxflux = matchxflux[alreadymatched]
            alreadymatchedrflux = matchrflux[alreadymatched]

            alreadymatcheddist = match_7_dists[alreadymatched]
            if fullfl <= alreadymatchedxflux:
                continue
            elif fullfl > alreadymatchedxflux:
                match_7.remove(good_gals)
                matchxflux.remove(alreadymatchedxflux)
                matchrflux.remove(alreadymatchedrflux)
                match_7_dists.remove(alreadymatcheddist)
                matchedind.remove(matchedind[alreadymatched])
                matchdec_diff.remove(matchdec_diff[alreadymatched])
                matchra_diff.remove(matchra_diff[alreadymatched])

        match_7.append(good_gals)
        match_7_dists.append(currentdist)
        matchxflux.append(fullfl)
        matchedind.append(i)
        matchrflux.append(galrflux)
        matchra_diff.append(delta_ra_good)
        matchdec_diff.append(delta_dec_good)
    return np.array(match_7), np.array(match_7_dists), np.array(matchxflux), np.array(matchedind), np.array(matchrflux), np.array(matchra_diff), np.array(matchdec_diff)

def magdifftointenrat(magdiff):
    '''
    converts magnitude difference to intensity ratio
    '''
    return 100**(magdiff/5)

def catmatch(skycoord, catcoords):
    '''
    uses astropy skycoord, just takes nearest position
    not advised for use
    '''
    id1, d2d1, d3d1 = skycoord.match_to_catalog_sky(catcoords)
    #id1 correspond to the GSW indices
    #id1[good1] correspond to the good gsw indices
    good1 = np.where(np.array(d2d1/u.degree) < 7*arcsec)[0] #these will correspond to 3xmm indices
    return id1, good1, d2d1
os.chdir(r'C:\Users\mugdhapolimera\github\xray')
xmm3 = pd.read_csv('3XMM_DR8_slim.csv')
ra1,dec1 = np.array(xmm3['sc_ra']), np.array(xmm3['sc_dec'])
fullflux = 0.9*(xmm3.sc_ep_2_flux + xmm3.sc_ep_3_flux + xmm3.sc_ep_4_flux + 
                xmm3.sc_ep_5_flux)
full = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\ECO+RESOLVE_filter_new.csv')
ra2,dec2 = np.array(full['radeg']), np.array(full['dedeg'])
rfluxes = np.array(full['absrmagcorr'])

#ra2, dec2 = np.loadtxt(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\RESOLVE_filter_coords2.txt',unpack = True) 
#ra1, dec1 = np.genfromtxt(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\xray\xmm3dr8.txt',usecols=(3,4),skip_header = 1,unpack = True) 
goodm1 = len(ra1)
match_7,match_7_dists,matchxflux,matchedind,matchrflux,matchra_diff,matchdec_diff =catmatch_act(ra1, dec1, ra2, dec2, rfluxes, fullflux, goodm1)
print (match_7)
full_xraymatched = full.iloc[match_7]
xmm_xraymatched = xmm3.iloc[matchedind]

#Convert X-ray flux (erg/cm^2/ Mpc) to Luminosity 
# L = 4*pi*d^2 * fx
# d = cz/H0
xrayagnprop = pd.read_csv('xmmagnprop.csv')

H0 = 70 #km/s/Mpc
Mpc = 3.086e24 #cm
d = (full_xraymatched.cz/H0) * Mpc #cm

Lx = np.log10(4*3.14*np.array(d**2) * np.array(fullflux[matchedind]))
plt.plot(Lx,np.log10(full_xraymatched.sfr_nuv), 'o')
xaxis = np.linspace(min(Lx)-1, max(Lx)+1)
plt.plot(xaxis, xaxis-40+np.log10(0.66),'b-.')
plt.plot(xaxis, xaxis-40+np.log10(0.66)-0.6,'g-.')
plt.xlim(38,44)
plt.ylim(-2,4)
Lx_to_sfr = Lx-40+np.log10(0.66)-0.6
xrayagn = np.log10(full_xraymatched.sfr_nuv) < Lx_to_sfr

full_xraymatched['logxray_lum'] = Lx
full_xraymatched['xmm_id'] = np.array(xmm_xraymatched['srcid'])
full_xraymatched['xrayagn'] = xrayagn
xmm_xraymatched.to_csv('XMM3_matched.csv')
full_xraymatched.to_csv('ECO+RESOLVE_xray_new.csv')

flags = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\eco+resolve_emlineclass_filter.csv')
flags.index = flags.galname
flags_xray = flags.loc[full_xraymatched.name]
print(flags_xray.iloc[np.where(flags_xray['sftoagn'])[0]])

plt.figure()
plt.hist(np.array(xrayagnprop.lgm_tot_p50[xrayagnprop.lgm_tot_p50 > 0]),
                  bins = 'fd')
yaxis= np.arange(0,60)
plt.plot(9.5*np.ones(len(yaxis)), yaxis, 'k-.')