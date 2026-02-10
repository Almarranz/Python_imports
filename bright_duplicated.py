#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  7 15:05:28 2026

@author: amartinez

 Remove duplicate/ambiguous bright sources based on angular proximity.

 Logic:
   1) Select "bright" sources with table[mag_col] < mag_cut.
   2) Find all pairs among bright sources within 'threshold' on sky.
   3) Any source participating in one or more close pairs is flagged as duplicate.
   4) Keep only bright sources that are NOT in close pairs; keep all non-bright sources.
   5) Return filtered table, mask (True=kept), and (optionally) removed indices.

 Parameters
 ----------
 table : array-like
     Catalog with magnitude column `mag_col` and length equal to `coords`.
     Typically an `astropy.table.Table` or `pandas.DataFrame`.
 coords : astropy.coordinates.SkyCoord
     Sky coordinates, length must match `len(table)`.
 band : str,
     Name of the magnitude column in `table`, default 'H'.
 mag_lim : float
     Magnitude threshold; sources with mag < mag_cut are considered "bright".
 sep_lim : astropy.units.Quantity
     Angular separation threshold to flag duplicates (default 0.5 arcsec).
 
     
 """


def bright_dup(table, coord, band, mag_lim, sep_lim):
    
    import numpy as np
    from astropy import units as u
    from astropy.coordinates import search_around_sky

    
    
    
# 1. Select bright stars
    if table[band].unit.is_equivalent(u.mag):
        bright_mask = table[band].value < mag_lim    
    else:
        bright_mask = table[band] < mag_lim
    bright_coords = coord[bright_mask]
    
    # 2. Search for all pairs of bright stars within a distance
    threshold = sep_lim * u.arcsec
    
    idx1, idx2, sep2d, _ = search_around_sky(
            bright_coords,
            bright_coords,
            threshold
    )
    
    # check = np.c_[idx1,idx2,sep2d.to(u.mas).value]
    
    # 3. Remove self-matches (each star matches itself at sep=0)
    self_matches = idx1 != idx2
    idx1 = idx1[self_matches]
    idx2 = idx2[self_matches]
    
    # 4. Any star involved in a "close pair" is considered a duplicate
    duplicate_bright = np.zeros(len(bright_coords), dtype=bool)
    duplicate_bright[np.unique(idx1)] = True
    duplicate_bright[np.unique(idx2)] = True
    
    # 5. We keep only the isolated bright stars
    keep_bright = np.logical_not(duplicate_bright)
    
    # 6. Build final mask for the full catalog
    final_mask = np.ones(len(table), dtype=bool)
    final_mask[bright_mask] = keep_bright
    
    # 7. Apply filtering
    table = table[final_mask]
    
    return table

# from typing import Tuple, Sequence, Optional
# import numpy as np
# import astropy.units as u
# from astropy.coordinates import SkyCoord
# from astropy.table import Table
# from astroquery.utils import commons  # only for unit hygiene (optional)
# from astropy.coordinates.search_around_sky import search_around_sky

# def remove_close_duplicates(
#     table,
#     coords: SkyCoord,
#     mag_col: str = 'H',
#     mag_cut: float = 13.0,
#     threshold: u.Quantity = 0.5 * u.arcsec,
#     return_indices: bool = True
# ) -> Tuple:
#     """
#     Remove duplicate/ambiguous bright sources based on angular proximity.

#     Logic:
#       1) Select "bright" sources with table[mag_col] < mag_cut.
#       2) Find all pairs among bright sources within 'threshold' on sky.
#       3) Any source participating in one or more close pairs is flagged as duplicate.
#       4) Keep only bright sources that are NOT in close pairs; keep all non-bright sources.
#       5) Return filtered table, mask (True=kept), and (optionally) removed indices.

#     Parameters
#     ----------
#     table : array-like
#         Catalog with magnitude column `mag_col` and length equal to `coords`.
#         Typically an `astropy.table.Table` or `pandas.DataFrame`.
#     coords : astropy.coordinates.SkyCoord
#         Sky coordinates, length must match `len(table)`.
#     mag_col : str, optional
#         Name of the magnitude column in `table`, default 'H'.
#     mag_cut : float, optional
#         Magnitude threshold; sources with mag < mag_cut are considered "bright".
#     threshold : astropy.units.Quantity, optional
#         Angular separation threshold to flag duplicates (default 0.5 arcsec).
#     return_indices : bool, optional
#         If True, also return indices of removed rows (relative to input `table`).

#     Returns
#     -------
#     filtered_table
#         Same type as `table` (if Table/DataFrame, slicing preserves type).
#     mask : np.ndarray (bool)
#         Boolean mask of length len(table): True for kept entries, False for removed.
#     removed_idx : np.ndarray (int), optional
#         Indices (0-based) of rows removed; returned only if `return_indices=True`.

#     Notes
#     -----
#     - Self-matches are excluded.
#     - This treats *any* member of a close pair as a duplicate (conservative).
#       If you prefer “one representative per pair”, you could instead pick, e.g.,
#       the brightest or highest-S/N within each close group.

#     Examples
#     --------
#     >>> # Given gns1 (Table) and gns1_c (SkyCoord):
#     >>> filtered, mask, removed = remove_close_duplicates(gns1, gns1_c, 'H', 13, 0.5*u.arcsec)
#     >>> # Save or continue with `filtered`
#     """
#     # Basic checks
#     n = len(table)
#     if len(coords) != n:
#         raise ValueError(f"`coords` length ({len(coords)}) must match `table` length ({n}).")
#     if mag_col not in table.colnames if isinstance(table, Table) else (mag_col in table.columns):
#         raise KeyError(f"Column '{mag_col}' not found in table.")

#     # 1) Bright selection
#     mags = table[mag_col]
#     bright_mask = (mags < mag_cut)
#     n_bright = np.count_nonzero(bright_mask)
#     if n_bright == 0:
#         # Nothing to filter; return original
#         mask = np.ones(n, dtype=bool)
#         return (table, mask, np.array([], dtype=int)) if return_indices else (table, mask)

#     bright_coords = coords[bright_mask]

#     # 2) Find pairs among bright sources within threshold
#     #    `search_around_sky` returns indices into the *input arrays* you pass in.
#     #    Here both inputs are bright_coords, so idx1/idx2 are 0..(n_bright-1).
#     idx1, idx2, sep2d, _ = search_around_sky(
#         bright_coords, bright_coords, threshold
#     )

#     # 3) Remove self-matches (i.e., i==j at sep=0)
#     not_self = (idx1 != idx2)
#     idx1 = idx1[not_self]
#     idx2 = idx2[not_self]
#     # (Optional) also keep only one of (i,j) / (j,i) by enforcing idx1 < idx2
#     # This doesn't change which stars are flagged, but avoids double counting
#     order = idx1 < idx2
#     idx1 = idx1[order]
#     idx2 = idx2[order]

#     # 4) Any source involved in a close pair is a duplicate (conservative)
#     duplicate_bright = np.zeros(n_bright, dtype=bool)
#     if idx1.size > 0:
#         duplicate_bright[np.unique(idx1)] = True
#         duplicate_bright[np.unique(idx2)] = True

#     # 5) Keep only isolated bright sources
#     keep_bright = ~duplicate_bright

#     # 6) Build final mask for the full table
#     final_mask = np.ones(n, dtype=bool)
#     final_mask[bright_mask] = keep_bright

#     # 7) Apply filtering
#     # Preserve input type on slicing
#     filtered = table[final_mask]

#     if return_indices:
#         removed_idx = np.flatnonzero(~final_mask)
#         return filtered, final_mask, removed_idx
#     else:
#         return filtered, final_mask