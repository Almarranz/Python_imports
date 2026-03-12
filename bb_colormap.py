#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 22:13:44 2026

@author: amartinez
"""
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import LogNorm
import numpy as np
from astromy_ds9 import ds9_norm   # assuming this is your local module
import IPython
import re
from astropy import units as u

def load_ds9_pseudocolor_lut(path, n=256):
    """
    Parse a DS9 'PSEUDOCOLOR' LUT file that defines per-channel control points:
        RED:   (x0,y0)(x1,y1)...
        GREEN: (x0,y0)(x1,y1)...
        BLUE:  (x0,y0)(x1,y1)...

    Returns
    -------
    cmap : matplotlib.colors.ListedColormap
        A Matplotlib colormap approximating DS9's linear interpolation between control points.
    """
    with open(path, 'r', encoding='utf-8') as f:
        txt = f.read()

    # Extract control points for a given channel name
    def extract(channel):
        # e.g. "RED:\n(0,0)(0.5,1)(1,1)"
        m = re.search(rf'{channel}:\s*([()\d\.\,\s]+)', txt, flags=re.IGNORECASE)
        if not m:
            raise ValueError(f'Channel "{channel}" not found in LUT file.')
        pairs = re.findall(r'\(([^)]+)\)', m.group(1))  # ['0,0', '0.5,1', '1,1', ...]
        arr = np.array([[float(a), float(b)] for a, b in (p.split(',') for p in pairs)], dtype=float)
        # Sort by x (position)
        arr = arr[np.argsort(arr[:, 0])]
        return arr

    red_cp   = extract('RED')
    green_cp = extract('GREEN')
    blue_cp  = extract('BLUE')

    # If values appear to be in 0..255, normalize to 0..1
    for arr in (red_cp, green_cp, blue_cp):
        if arr[:, 0].max() > 1.0 or arr[:, 1].max() > 1.0:
            arr[:, 0] = arr[:, 0] / arr[:, 0].max()
            arr[:, 1] = arr[:, 1] / 255.0

    x = np.linspace(0.0, 1.0, n)
    r = np.interp(x, red_cp[:,   0], red_cp[:,   1])
    g = np.interp(x, green_cp[:, 0], green_cp[:, 1])
    b = np.interp(x, blue_cp[:,  0], blue_cp[:,  1])

    lut = np.stack([r, g, b], axis=1)
    return ListedColormap(lut, name='ds9_bb')

# =============================================================================
# HOW to USE IT
# Build the DS9 'bb' colormap from your LUT file
# =============================================================================

# bb_lut_path = "/Users/amartinez/Desktop/PhD/images/colormaps/bb.lut" # This can be donwloaded from ds9. Color, colomap param, file, save
# cmap_bb = load_ds9_pseudocolor_lut(bb_lut_path, n=256)