#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 10:48:15 2026

@author: amartinez
"""

from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import numpy as np

def l_function(table, band, b_with, m_l = None, m_h = None, min_counts = None):
    
    if m_l == None:
        m_l = np.min(table[band].value)
        m_h = np.max(table[band].value)
        
    num_bins = np.arange(m_l, m_h ,b_with )
    # Create the histogram
    mag = table[band].value
    counts, bin_edges = np.histogram(mag, bins=num_bins)
    
    if min_counts is None:
        min_counts = 1
    c_mask = counts > min_counts

    bins = (bin_edges[:-1] + bin_edges[1:]) / 2

    counts = counts[c_mask]
    bins = bins[c_mask]
    
    return bins, counts


def plot_two_pm_hists(table, var1, var2, sym1, sym2, *,
                      bins=50, figsize=(8,4), title1 = None,  title2 = None):
    """
    Plot side-by-side histograms for two columns, with mean and std in the legend.

    Parameters
    ----------
    table : astropy.table.Table or pandas.DataFrame
        The catalog containing the data.
    var1, var2 : str
        Column names to plot (e.g., 'pmra', 'pmde').
    sym1, sym2 : str
        LaTeX math symbols WITHOUT dollar signs, e.g., r'\\mu_{RA}', r'\\mu_{Dec}'.
        These are used in axis labels and in the overline for the mean.
    bins : int, sequence, str or tuple, optional (default: 'auto')
        Bin specification for the histograms.
        - If a single value (e.g., 50 or 'auto'), it is applied to both histograms.
        - If a 2-tuple (bins1, bins2), each is applied to the corresponding panel.
    figsize : 2-tuple, optional (default: (10, 4))
        Figure size in inches, e.g., (width, height).

    Notes
    -----
    - Pass symbols like r'\\mu_{RA}', not '$...$'. The function adds the dollar signs.
    - Legends show \\overline{symbol} (mean) and \\sigma (std) with 2 decimals.
    """
    # Normalize bins argument to a pair (bins1, bins2)
    if isinstance(bins, tuple) and len(bins) == 2:
        bins1, bins2 = bins
    else:
        bins1 = bins2 = bins

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    if title1 is not None:
        ax.set_title(title1)
    if title2 is not None:
        ax2.set_title(title2)

    # ---- FIRST ----
    vals1 = np.asarray(table[var1])
    ax.hist(vals1, bins=bins1, histtype='step',
            label=fr'$\overline{{{sym1}}} = {np.mean(vals1):.2f}$' +'\n'
                  + fr'$\sigma = {np.std(vals1):.2f}$')
    ax.set_xlabel(fr'${sym1}$ [mas/yr]')
    ax.set_ylabel('# stars')
    ax.legend()

    # ---- SECOND ----
    vals2 = np.asarray(table[var2])
    ax2.hist(vals2, bins=bins2, histtype='step',
             label=fr'$\overline{{{sym2}}} = {np.mean(vals2):.2f}$' +'\n'
                   + fr'$\sigma = {np.std(vals2):.2f}$')
    ax2.set_xlabel(fr'${sym2}$ [mas/yr]')

    ax2.legend()
    plt.tight_layout()
    plt.show()
    
def sig_cl(x, y,s):
    mx, lx, hx = sigma_clip(x , sigma = s, masked = True, return_bounds= True)
    my, ly, hy = sigma_clip(y , sigma = s, masked = True, return_bounds= True)
    m_xy = np.logical_and(np.logical_not(mx.mask),np.logical_not(my.mask))
    
    return m_xy, [lx,hx,ly,hy]