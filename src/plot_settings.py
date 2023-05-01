#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 17:11:00 2023

@author: jaya
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pickle
# from floquet_qubit import Floquet1dQubit
from wcs_colors import ccycle_floquet as ccycle
import os
import scipy
import scipy.constants
import seaborn as sns
from collections import OrderedDict
import matplotlib.pyplot as plt
#%%

##### setting figure size #####

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

#rc('text', usetex=True)
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'


#%%
##### initializing colors #####

'''
Physical constants
'''
pi = np.pi
h = scipy.constants.h
hbar = scipy.constants.hbar
e = scipy.constants.elementary_charge
phi0 = h/(2*e)
reducedphi0 = hbar/(2*e)
'''
Initialize Colors
'''
zero = sns.color_palette("RdBu", 10)[0]
one = sns.color_palette("RdBu", 10)[-1]
two = 'green'
two = (93./255, 136./255, 54./255)
three = (141/255., 52/255., 120/255.)
four = sns.color_palette("Paired", 10)[7]
four = 'darkorange'
six = (236/255., 225/255., 51/255.)
five =    (24./255., 157./255., 194./255)
seven = 'palevioletred'#(220/255., 148/255., 162/255.)
eight = (169/255., 198/255., 105/255.)
nine = (166/255., 156/255., 199/255.)
ten = sns.color_palette("YlOrRd", 10)[2]
eleven = (142/255., 194./255, 198./255) 
twelve = 'plum'#sns.color_palette("PiYG", 10)[3]
thirteen = sns.color_palette("RdYlGn", 10)[-5]
fourteen = sns.color_palette("RdYlBu", 10)[-4]
fifteen = sns.color_palette("RdPu", 10)[2]
sixteen = sns.light_palette("navy", reverse=True)[-1]
arbit_high= sns.color_palette("Set2")[-2]
arbit_higher = 'lightgray'
'''
Potentials
'''
LeatherJacket = '#708090'
BlackBean = '#32174D'#(61/255., 12/255., 2/255.)

states = [
        zero,
        one,
        two,
        three,
        four,
        five,
        six,
        seven,
        eight,
        nine,
        ten,
        eleven,
        twelve,
        thirteen,
        fourteen,
        fifteen,
        sixteen,
        arbit_high,
        arbit_higher
        ]
states_all = states + [states[-1]]*100 
#%%
##### plot parameters #####

if 1: # plot parameters PRB template
    mpl.rc('font',family='serif')
#    mpl.rcParams['font.serif'] = 'Times New Roman'
    mpl.rc('font',size=9)
    mpl.rc('font',size=9)
    mpl.rcParams["text.usetex"]=True
    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams['font.family'] = 'STIXGeneral'
    mpl.rcParams['axes.formatter.useoffset'] = False
    
    mpl.rcParams['ytick.right'] = False
    mpl.rcParams['xtick.top'] = False
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['legend.fontsize'] = 9
    mpl.rcParams['legend.frameon'] = False
    mpl.rcParams['axes.labelsize'] = 9
    mpl.rcParams['xtick.labelsize'] = 9
    mpl.rcParams['ytick.labelsize'] = 9
    mpl.rcParams['xtick.major.pad']=  2 #3.5
    mpl.rcParams['ytick.major.pad']=  2 #3.5
    mpl.rcParams['axes.labelpad'] = 1 #4.0
    mpl.rcParams['legend.handlelength'] = 1.0#2.0
    mpl.rcParams['legend.handletextpad'] = 0.4# 0.8
    mpl.rcParams['legend.columnspacing'] = 1.2# 2.0,
    mpl.rcParams['lines.markersize'] = 4.0
    mpl.rcParams["figure.figsize"] = [7.4, 3.5]