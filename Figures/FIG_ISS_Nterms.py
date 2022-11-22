#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Determine the number of terms included in the evaluation of the infinite-series
solution (ISS) of the Mason-Weaver equation.
"""

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('..') # make masonweaver_analytic available here

import masonweaver_analytic
from masonweaver_analytic import MW_adim_ISS

DRAFT = True

#%% CALCULATE

if DRAFT:
    Ntau = 100 # quick generate
else:
    Ntau = 1000 # hi-res
Ny = 100

alphas = [0.02, 0.1, 0.5, 2.0, 10., 100.]
Nterms = np.zeros((Ntau, len(alphas)), dtype = np.uint)
tau = np.logspace(-3.0, 3.0, Ntau)
y = np.linspace(0, 1, Ny)

for iy, alpha in enumerate(alphas):
    for ix in range(Ntau):
        # evaluate the concentration profile ISS of the MWE at tau, alpha
        cc = MW_adim_ISS(y, tau[ix], alpha)
        # retrieve the number of terms used in its calculation
        Nterms[ix, iy] = masonweaver_analytic.ISS_Nterms_latest
        

        
#%% PLOT
        
from matplotlib import rcParams
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['axes.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['figure.figsize'] = (3.149, 2.362) # 80mm x 50mm
rcParams['figure.subplot.left'] = 0.19  
rcParams['figure.subplot.right'] = 0.94    
rcParams['figure.subplot.bottom'] = 0.19
rcParams['figure.subplot.top'] = 0.94
rcParams['legend.frameon'] = False

fig = plt.figure(1)
plt.clf()
for iy, alpha in enumerate(alphas):
    if alpha < 10:
        lbl = '$\\alpha = {0:.2f}$'.format(alpha)
    else:
        lbl = '$\\alpha = {0:.0f}$'.format(alpha)
    plt.loglog(tau, Nterms[:,iy], label = lbl)
fig.axes[0].set_xticks([0.001, 0.01, 0.1, 1.0, 10, 100, 1000])
plt.legend(labelspacing =  0.2)
plt.ylabel('$N_{\\rm terms}$ (ISS)')
plt.xlabel('$\\tau$')

if not DRAFT:
    plt.savefig('FIG_ISS_Nterms.png', dpi = 450)

plt.show()
