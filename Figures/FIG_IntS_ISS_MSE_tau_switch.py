#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Caclulate MSE between IntS and ISS, plot tau_switch
"""

import numpy as np
from numpy import log10

import matplotlib.pyplot as plt

import sys
sys.path.append('..') # make masonweaver_analytic available here

from masonweaver_analytic import MW_adim_ISS, MW_adim_IntS
from masonweaver_analytic import MSE, tau_switch

#%% SET UP

DRAFT = True
SAVEFIG = False
outfn = 'FIG_IntS_ISS_MSE_tau_switch.png'

# set resolution
if DRAFT:
    Nalphas = 30
    Ntaus = 25
    Ny = 100
else:
    Nalphas = 120
    Ntaus = 120
    Ny = 100


#%% CALCULATE

# scan a grid, calculate MSE
y = np.linspace(0, 1, Ny)
alphas = np.logspace(log10(0.02), log10(20.0), Nalphas)
taus = np.logspace(-3, 0, Ntaus)
eps = np.zeros((Ntaus, Nalphas))

for iy, alpha in enumerate(alphas):
    for ix, tau in enumerate(taus):
        cISS = MW_adim_ISS(y, tau, alpha)
        cIntS = MW_adim_IntS(y, tau, alpha)
        eps[ix, iy] = MSE(cIntS, cISS)


#%% PLOT

from matplotlib import rcParams
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['axes.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['figure.figsize'] = (3.149, 3.149*0.8) # 80mm x xxx mm
rcParams['figure.subplot.left'] = 0.19  
rcParams['figure.subplot.right'] = 0.9   
rcParams['figure.subplot.bottom'] = 0.19
rcParams['figure.subplot.top'] = 0.94
rcParams['legend.frameon'] = False

plt.figure(1)
plt.clf()
plt.pcolormesh(log10(alphas), log10(taus), log10(eps))
plt.colorbar(label='$^{10}\\log({\\rm MSE})$')
plt.xlabel('$^{10}\\log (\\alpha)$')
plt.ylabel('$^{10}\\log (\\tau)$')

aa = np.logspace(log10(0.02), log10(20.0), 1000)
tswitch = [tau_switch(a) for a in aa]
plt.plot(log10(aa), log10(tswitch), 'r')
#plt.plot([log10(aa[-1]), log10(aa[-1])], [log10(tswitch[-1]), -3.0], 'r')

if SAVEFIG:
    plt.savefig(outfn, dpi = 450)
    
plt.show()

