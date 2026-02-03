#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 22:19:05 2025

@author: amartinez
"""


from matplotlib import rc
from matplotlib import rcParams
import matplotlib.pyplot as plt
import IPython
def set_style():
    rcParams.update({
        "text.usetex": False,
        "font.family": "sans",
        "font.sans-serif": ["Palatino"]
    })

    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    rc("font", **{"family": "serif", "serif": ["Palatino"]})
    # %%plotting parametres
    rcParams.update({'xtick.major.pad': '7.0'})
    rcParams.update({'xtick.major.size': '7.5'})
    rcParams.update({'xtick.major.width': '1.5'})
    rcParams.update({'xtick.minor.pad': '7.0'})
    rcParams.update({'xtick.minor.size': '3.5'})
    rcParams.update({'xtick.minor.width': '1.0'})
    rcParams.update({'ytick.major.pad': '7.0'})
    rcParams.update({'ytick.major.size': '7.5'})
    rcParams.update({'ytick.major.width': '1.5'})
    rcParams.update({'ytick.minor.pad': '7.0'})
    rcParams.update({'ytick.minor.size': '3.5'})
    rcParams.update({'ytick.minor.width': '1.0'})
    rcParams.update({'font.size': 20})
    rcParams.update({'figure.figsize':(10,5)})
    rcParams.update({
        "text.usetex": False,
        "font.family": "sans",
        "font.sans-serif": ["Palatino"]})
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'
    rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rcParams.update({'figure.max_open_warning': 0})# a warniing for matplot lib pop up because so many plots, this turining it of
    # Enable automatic plotting mode
    # IPython.get_ipython().run_line_magic('matplotlib', 'auto')
    IPython.get_ipython().run_line_magic('matplotlib', 'inline')