#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:52:39 2026

@author: amartinez
"""





# =============================================================================
# MAIN FUNCTION (from: https://spacetelescope.github.io/hst_notebooks/notebooks/WFC3/zeropoints/zeropoints.html)
# =============================================================================
def get_vegazp(mjd, band, aper=6.0):
    
    import os
    import tarfile
    import numpy as np
    from astropy.time import Time
    
    # =============================================================================
    # CONFIGURATION — set once and reuse
    # =============================================================================
    CDBS_TAR = "/Users/amartinez/Desktop/Projects/SOMA_HST_pm/synphot_cdbs/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar"
    CDBS_DIR = "/Users/amartinez/Desktop/Projects/SOMA_HST_pm/synphot_cdbs/"
    TRDS_DIR = CDBS_DIR + "hlsp_reference-atlases_hst_multi_everything_multi_v11_sed/grp/redcat/trds/"

    
    
    from synphot import Observation


    
    
    # Make sure CDBS is extracted only once
    def _ensure_cdbs():
        if not os.path.exists(TRDS_DIR):
            print("Extracting CDBS throughput library…")
            with tarfile.open(CDBS_TAR, "r") as tar:
                tar.extractall(path=CDBS_DIR, filter="data")
        os.environ["PYSYN_CDBS"] = TRDS_DIR
    
    """
    Compute the VEGAMAG zeropoint for HST/WFC3-IR for a given date and filter.

    Parameters
    ----------
    mjd : float or str
        Observation date in MJD or ISO (converted automatically)
    band : str
        HST WFC3/IR filter name (e.g. 'f160w', 'f125w', 'f128n')
    aper : float
        Photometric aperture radius in arcsec (6.0 = infinite)

    Returns
    -------
    ZP_Vega : float
        VEGAMAG zeropoint for the specified filter and date
    """

    # --- Ensure CDBS is available ---
    _ensure_cdbs()

    # --- Normalize date ---
    try:
        mjd = float(mjd)
    except:
        mjd = Time(mjd).mjd
    
    import stsynphot as stsyn
    # --- Set Vega spectrum (once per session) ---
    vega_url = "https://ssb.stsci.edu/trds/calspec/alpha_lyr_stis_010.fits"
    if not hasattr(stsyn, "Vega") or stsyn.Vega is None:
        stsyn.Vega = stsyn.spectrum.SourceSpectrum.from_file(vega_url)

    # --- Build obsmode ---
    obsmode = f"wfc3,ir,{band},mjd#{mjd},aper#{aper}"

    # --- Bandpass ---
    bp = stsyn.band(obsmode)

    # --- PHOTFLAM (inverse sensitivity) ---
    photflam = bp.unit_response(stsyn.conf.area)

    # --- STMAG zeropoint ---
    stmag = -21.1 - 2.5 * np.log10(photflam.value)

    # --- pivot wavelength ---
    photplam = bp.pivot()

    # --- ABMAG zeropoint ---
    abmag = stmag - 5 * np.log10(photplam.value) + 18.6921

    # --- VEGAMAG zeropoint: Vega has mag=0 in VEGAMAG ---
    obs = Observation(stsyn.Vega, bp, binset=bp.binset)
    ZP_Vega = obs.effstim(flux_unit="obmag", area=stsyn.conf.area)

    return ZP_Vega.value*-1