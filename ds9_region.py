#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
def region(table, coor1, coor2,
                     name='regions',
                     save_in='.',
                     color='green',
                     wcs='fk5',
                     marker='cross'):
"""


def region(table, coor1, coor2,
                     name='regions',
                     save_in='.',
                     color='green',
                     wcs='fk5',
                     marker='cross'):
    
    import os
    import astropy.units as u

    """
    Writes a DS9 region file from a catalog containing RA/Dec coordinates.

    Parameters
    ----------
    table : astropy Table or pandas DataFrame
        Catalog with RA/Dec columns.
    ra_col, dec_col : str
        Column names for RA and Dec (in degrees).
    name : str
        Output filename (without extension).
    save_in : str
        Directory where the .reg file will be saved.
    color : str or None
        DS9 color name. Default: 'green'.
    wcs : str or None
        WCS system to write (e.g. 'fk5', 'icrs', 'galactic').
        If None → 'physical'.
    marker : str
        DS9 marker style (e.g. 'cross', 'circle', 'box').
    """
    

    # Determine final color
    final_color = color if color is not None else 'green'

    # File path
    filepath = os.path.join(save_in, f'{name}.reg')

    with open(filepath, 'w') as regfile:
        # Header
        regfile.write('# Region file format: DS9 version 4.1\n')
        regfile.write(
            f'global color={final_color} dashlist=8 3 width=1 '
            'font="helvetica 10 normal roman" select=1 highlite=1 '
            'dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        )

        # WCS / physical system
        if wcs is None:
            regfile.write('physical\n')
        else:
            regfile.write(f'{wcs}\n')
        
        unidades = table[coor1].unit
        # Write all points
        if unidades is not None:
        
            if marker == 'circulos':        
                for ra, dec in zip(table[coor1].value, table[coor2].value):
                    regfile.write(f'circle({ra},{dec}, 100")\n')
                for ra, dec in zip(table[coor1].value, table[coor2].value):
                    regfile.write(f'point({ra},{dec}) # point= cross\n')
                    
            else:
                for ra, dec in zip(table[coor1].value, table[coor2].value):
                    regfile.write(f'point({ra},{dec}) # point={marker}\n')
        else: 
            if marker == 'circulos':        
                for ra, dec in zip(table[coor1], table[coor2]):
                    regfile.write(f'circle({ra},{dec}, 100")\n')
                for ra, dec in zip(table[coor1], table[coor2]):
                    regfile.write(f'point({ra},{dec}) # point= cross\n')
                    
            else:
                for ra, dec in zip(table[coor1], table[coor2]):
                    regfile.write(f'point({ra},{dec}) # point={marker}\n')

    print(f"DS9 region file saved to: {filepath}")

def region_vectors(table,
                   ra_col, dec_col,
                   pmra_col, pmdec_col,
                   name='vector_regions',
                   save_in='.',
                   color='green',
                   wcs='fk5',
                   scale=10.0,
                   width = 1):
    """
    Write a DS9 region file containing proper-motion vectors.

    Parameters
    ----------
    table : astropy Table
        Table containing RA/Dec and proper motion components.
    ra_col, dec_col : str
        Column names for RA, Dec in degrees.
    pmra_col, pmdec_col : str
        Column names for proper motion components (e.g. mas/yr).
    name : str
        Output filename (without .reg extension).
    save_in : str
        Directory where the .reg is saved.
    color : str
        DS9 color.
    wcs : str or None
        WCS header string ('fk5', 'icrs', etc.), None → 'physical'.
    scale : float
        A scaling factor for vector lengths (in arcsec). Default = 10.
        Length is:  scale * sqrt(pmRA^2 + pmDec^2)
    """

    import os
    import numpy as np

    # Build output path
    filepath = os.path.join(save_in, f"{name}.reg")

    with open(filepath, 'w') as regfile:
        # Header
        regfile.write('# Region file format: DS9 version 4.1\n')
        regfile.write(
            f'global color={color} dashlist=8 3 width=1 '
            'font="helvetica 10 normal roman" select=1 highlite=1 '
            'dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        )

        # Coordinate system
        if wcs is None:
            regfile.write('physical\n')
        else:
            regfile.write(f'{wcs}\n')

        # Extract columns
        RA   = table[ra_col]
        DEC  = table[dec_col]
        PMRA = table[pmra_col]
        PMDEC = table[pmdec_col]

        # Handle units if present
        if hasattr(RA, 'unit'):   RA = RA.value
        if hasattr(DEC, 'unit'):  DEC = DEC.value
        if hasattr(PMRA, 'unit'): PMRA = PMRA.value
        if hasattr(PMDEC,'unit'): PMDEC = PMDEC.value

        # Compute vector lengths (in arcsec)
        # You may change to arcmin or degrees if preferred
        length = scale * np.sqrt(PMRA**2 + PMDEC**2)

        # Compute angle using atan2
        # atan2(y, x) = atan(PMDEC/PMRA) automatically handles quadrants
        angle_deg = np.degrees(np.arctan2(PMDEC, PMRA))

        # Write each vector
        for ra, dec, L, ang in zip(RA, DEC, length, angle_deg):
            regfile.write(f'# vector({ra},{dec},{L}\",{ang}) vector=1 width={width}\n')

    print(f"DS9 vector region file saved to: {filepath}")