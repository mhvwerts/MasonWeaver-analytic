#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test mass conservation of numerically evaluated analytic solutions to the 
Mason-Weaver equation.

This program scans a grid of (alpha, tau) values, selects the appropriate
non-dimensional solution function as laid out in the paper, and then
numerically integrates that analytic function using scipy.integrate.quad
(QUADPACK library [1][2]). The integral of the non-dimensional solution should
consistently evaluate as 1 if total mass of the system is conserved.


[1] https://doi.org/10.1007/978-3-642-61786-7
[2] https://netlib.org/quadpack/changes



Typical run times for this script observed on development system
    Intel Core i7 (12th gen., 2.4 GHz), Windows 10pro, 
    Python 3.9.13, numpy 1.22.3, scipy 1.9.1

DRAFT=True  noIntS=False          2.9 s
DRAFT=True  noIntS=True         140.5 s
DRAFT=False noIntS=False         42.9 s
DRAFT=False noIntS=True      not tested
"""

#%%

from time import process_time

import numpy as np
from numpy import log10

from scipy.integrate import quad

import matplotlib.pyplot as plt

import sys
sys.path.append('..') # masonweaver_analytic.py is in parent directory
import masonweaver_analytic
from masonweaver_analytic import MW_adim_ISS, MW_adim_IntS
from masonweaver_analytic import MW_adim_hiPe, MW_adim_hiPe2
from masonweaver_analytic import tau_switch


#%% SET UP

# If `noIntS` is set to True, the program will exclusively use ISS for
# alpha > 0.02 (alpha > masonweaver_analytic.ALPHA_SWITCH2) and will never use
# the faster IntS. This leads to a vast increase in computation time (50x)
# for this test, but will somewhat increase precision and mass conservation
# of the solution for a limited range of (alpha, tau).
noIntS = False 

# DRAFT enables a faster calculation on a less densely sampled grid
DRAFT = False


# set resolution
if DRAFT:
    Nalphas = 40
    Ntaus = 30
else:
    Nalphas = 120
    Ntaus = 150
    

#%% GET AND CHECK PARAMETERS

# Explicitly switch solutions
alpha_switch1 = masonweaver_analytic.ALPHA_SWITCH1 
alpha_switch2 = masonweaver_analytic.ALPHA_SWITCH2
alpha_switch3 = masonweaver_analytic.ALPHA_SWITCH3

# alphas should start above alpha_switch1, since hiPe2 does not conserve mass
alpha_start = 0.002
assert alpha_start > alpha_switch1,\
    'alphas should start above alpha_switch1, '\
    'since hiPe2 does not conserve mass'

#%% CALCULATE



# scan a grid, calculate integral over analytic solution
alphas = np.logspace(log10(0.002), 3, Nalphas)
taus = np.logspace(-3, 2, Ntaus)
totalmass = np.zeros((Ntaus, Nalphas))
funmap = np.zeros((Ntaus, Nalphas), dtype = np.uint)

t0 = process_time()

for iy, alpha in enumerate(alphas):
    for ix, tau in enumerate(taus):
        # select the appropriate function to be evaluated for mass conservation
        if alpha < alpha_switch1:
            fun = MW_adim_hiPe2;funcode=4
        elif alpha_switch1 <= alpha < alpha_switch2:
            fun = MW_adim_hiPe;funcode=3
        elif alpha_switch2 <= alpha < alpha_switch3:
            if noIntS: # use exclusively ISS
                fun = MW_adim_ISS;funcode=1
            else: # switch to faster IntS if appropriate
                tau_sw = tau_switch(alpha)
                if tau < tau_sw:
                    fun = MW_adim_IntS;funcode=2
                else:
                    fun = MW_adim_ISS;funcode=1
        else:
            fun = MW_adim_ISS;funcode=1
        funmap[ix, iy] = funcode
        # integrate the function using 'scipy.integrate.quad'
        # to obtain total mass for mass conservation test
        tot, abserr = quad(fun, 0, 1, args=(tau,alpha))
        totalmass[ix, iy] = tot

t1 = process_time()

print()
print('*'*70)
print('*** Mass conservation grid scan (tau, alpha) complete ***')
print('*'*70)

print('Highest value for mass integral: ', totalmass.max())
print('Lowest value for mass integral:  ', totalmass.min())

print('Grid calculation time:            {0:.2f} s ('.format(t1-t0)+
      ('pure ISS' if noIntS else 'ISS/IntS switch')+
      ')')

plt.figure(1)
plt.clf()
plt.pcolormesh(log10(alphas), log10(taus), totalmass)
plt.xlabel('$^{10}\log \alpha$')
plt.ylabel('$^{10}\log \tau$')
plt.colorbar(label = 'total mass integral')
plt.title('Mass conservation of analytic solution')

# draw a map of which function was used at the combination of (alpha, tau)
plt.figure(2)
plt.clf()
plt.pcolormesh(log10(alphas), log10(taus), funmap)
plt.xlabel('$^{10}\log \alpha$')
plt.ylabel('$^{10}\log \tau$')
plt.title('Solution map: 1=ISS, 2=IntS, 3=hiPe')
plt.colorbar(label = 'solution function')

plt.show()
