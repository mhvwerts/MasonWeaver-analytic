# -*- coding: utf-8 -*-
"""
Compare ISS, IntS and hiPe solutions at specific (tau, alpha) illustrating the
instability of ISS and IntS and the necessity of hiPe

alpha = 0.01
tau = 0.23
"""

SAVEFIG = False

#%% Imports

import numpy as np
import matplotlib.pyplot as plt


import sys
sys.path.append('..') # make masonweaver_analytic available here

from masonweaver_analytic import MW_adim_ISS, MW_adim_IntS, MW_adim_hiPe

#%% Configure plotting

from matplotlib import rcParams
rcParams['lines.linewidth'] = 1.0
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['axes.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['figure.figsize'] = (3.149, 3.149*1.0) # 80mm x xxx mm
rcParams['figure.subplot.left'] = 0.13
rcParams['figure.subplot.right'] = 0.95   
rcParams['figure.subplot.bottom'] = 0.13
rcParams['figure.subplot.top'] = 0.98
rcParams['legend.frameon'] = False


#%% Numerically evaluate different analytic solutions at same alpha, tau

alpha = 0.01
tau = 0.23

yy = np.linspace(0,1,1000)
ccISS = MW_adim_ISS(yy, tau, alpha)
ccInt = MW_adim_IntS(yy, tau, alpha)
cchiPe = MW_adim_hiPe(yy, tau, alpha)

#%% Plot results

plt.figure(1)
plt.clf()

plt.subplot(311)
plt.plot(yy,ccISS, 'r', label = 'ISS')
plt.gca().axes.xaxis.set_ticklabels([])
plt.legend(frameon = False)
plt.ylim(0,4)

plt.subplot(312)
plt.plot(yy,ccInt, 'b', label = 'IntS')
plt.gca().axes.xaxis.set_ticklabels([])
plt.ylabel('$c(y) / c_0$')
plt.legend(frameon = False)
plt.ylim(0,4)

plt.subplot(313)
plt.plot(yy,cchiPe, 'k', label = 'hiPe')
plt.legend(frameon = False)
plt.ylim(0,4)
plt.xlabel('$y$')

if SAVEFIG:
    plt.savefig('FIG_ISS_IntS_hiPe_instability.png',
                dpi = 350)

plt.show()
