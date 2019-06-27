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
            newdf.append([pcatalog.ra.value[place], pcatalog.dec.value[place], -99, -99, -99, pz[place],
                             -99, np.nan, np.nan, np.nan, np.nan,np.nan, np.nan, np.nan, np.nan, np.nan,np.nan,
                         np.nan])
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
                
            newdf.append([pcatalog.ra.value[place], pcatalog.dec.value[place], scatalog.ra.value[index],
                            scatalog.dec.value[index], d2d[index].arcsec, pz[place], sz[index],
                         sdss_mags["u"][index], sdss_mags["g"][index], sdss_mags["r"][index], 
                          sdss_mags["i"][index], sdss_mags["z"][index], bes_mags["U"][index], bes_mags["B"][index],
                          bes_mags["V"][index], bes_mags["R"][index], bes_mags["I"][index], s_names[index]])
    
    cnames = ["photRA", "photDec", "specRA", "specDec", "separation", "photZ", "specZ", "sdss_absmagu",
              "sdss_absmagg", "sdss_absmagr", "sdss_absmagi", "sdss_absmagz", 
              "bessell_absmagU", "bessell_absmagB", "bessell_absmagV", "bessell_absmagR", "bessell_absmagI",
             "specName"]
    return pd.DataFrame(newdf, columns = cnames)
