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

        # Write all points
        if table['l'].unit.is_equivalent(u.deg):
        
            if marker == 'circulos':        
                for ra, dec in zip(table[coor1].value, table[coor2].value):
                    regfile.write(f'circle({ra},{dec}, 100")\n')
                for ra, dec in zip(table[coor1].value, table[coor2].value):
                    regfile.write(f'point({ra},{dec}) # point= cross\n')
                    
            else:
                for ra, dec in zip(table[coor1], table[coor2]):
                    regfile.write(f'point({ra},{dec}) # point={marker}\n')
        else: 
            if marker == 'circulos':        
                for ra, dec in zip(table[coor1].value, table[coor2].value):
                    regfile.write(f'circle({ra},{dec}, 100")\n')
                for ra, dec in zip(table[coor1], table[coor2]):
                    regfile.write(f'point({ra},{dec}) # point= cross\n')
                    
            else:
                for ra, dec in zip(table[coor1], table[coor2]):
                    regfile.write(f'point({ra},{dec}) # point={marker}\n')

    print(f"DS9 region file saved to: {filepath}")


