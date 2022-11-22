#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example of the calculation of the analytic solution of the Mason-Weaver
equation.

This script calculates the same profiles as the `example.py` from the
`MasonWeaver-finite-diff` package, but uses analytic expressions instead of
a numerical finite-difference scheme.

L. Barthe and M. H. V. Werts, 2022
CNRS/ENS Rennes, lab. SATIE

published under the CeCILL v2.1 license
    https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
"""

import numpy as np
import matplotlib.pyplot as plt

from masonweaver_analytic import MW_c_profile

#
# Evaluate the analytic solution of the Mason-Weaver equation
# over 'Ntimesteps' time steps (exponentially increasing delta t)
#
# The parameters 'D' and 'sg' are for 60 nm diam. spherical gold nanoparticles
# in water at +4°C.
# Everything calculated using SI units: meters, seconds, ...
# For plotting, occasional conversions to centimeters, hours.
#

# Calculation parameters set according to `example.py` from
#   `MasonWeaver-finite-diff` to enable direct comparison

Ntimesteps = 500
Ngridpoints = 1000
t_end = 3.6e6

z_max = 0.01        # 1 cm cell height 
D = 4.32e-12        # for 60 nm gold nanosphere in water at 4 degC
sg = 2.47e-9*9.81   # for 60 nm gold nanosphere in water at 4 degC



# First, calculate the time points where the analytic solution will be 
# evaluated/plotted. We need to choose the same time points as in the
# example with the `MasonWeaver-finite-diff` package (exponentially expanding
# time grid).
# This calculation is done based on the code in that package.

t0 = D/(sg**2)
_tau_end = t_end/t0
_k_tau = np.log(_tau_end+1)/Ntimesteps
_tau = -1.0 + np.exp(_k_tau*np.arange(0, Ntimesteps+1))
tt = _tau*t0 # numpy vector containing the time points


###
# concentration profiles at specified time points

plt.figure(1)
plt.clf()
plt.title('Sedimentation 60nm diam. gold NP in water at 4°C')
for t in tt[::50]: # only every 50th time point
    solnz, solnc = MW_c_profile(t, z_max, D, sg, Nz = 500)
    plt.plot(solnz*100, solnc,
             label='t={:.0f} h '.format(t/3600))
plt.ylim(ymin = 0.0, ymax = 3.0)
plt.xlabel('height / cm')
plt.ylabel('rel. concentration')
plt.legend()


###
# time behaviour
#
#TODO: provide functionality in `mason_weaver_analytic` 
#      for calculating the concentration in a single
#      position in a sedimentation profile using non-adimensional units
#
# plt.figure(2)
# plt.clf()
# plt.title('Concentration at bottom of cell')
# z = 0.0 # the position of interest
# plt.plot(tt/3600, [TODO_MW_c(z, t, z_max, D, sg) for t in tt])
# plt.ylabel('rel. concentration')
# plt.xlabel('time / h')


###
# convergence to analytical steady-state

plt.figure(3)
plt.clf()
plt.title('Sedimentation equilibrium gradient')

# analytical expression steady-state from Midelet PPSC 2017
z0 = D/sg
B = z_max / (z0*(1-np.exp(-z_max/z0)))
z = np.linspace(0, z_max, 500)
plt.plot(z*100, B * np.exp(-z/z0),
         label = 'analytic expression')

solnz, solnc = MW_c_profile(tt[-1], z_max, D, sg, Nz = 500)
plt.plot(solnz*100, solnc,
         label = 'PWS very long time')

solnz, solnc = MW_c_profile(None, z_max, D, sg, Nz = 500)
plt.plot(solnz*100, solnc,
         label = 'PWS steady-state') # This is exactly the same function as
                                     # the analytic expression from Midelet
                                     # PPSC 2017

plt.xlabel('height / cm')
plt.ylabel('rel. concentration')
plt.legend()

plt.show()
