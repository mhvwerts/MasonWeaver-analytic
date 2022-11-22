#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot figure of the computation time for evaluation of IntS and ISS, together
with the mean-squared error (MSE) of the IntS with respect to the ISS.

The computation times are measured and the MSE calculated in the script
`FIG_IntS_ISS_comp_time_measure.py`. The present script reads the result from
a 'pickle' file and plots them.
"""

#%%
from matplotlib import pyplot as plt
import numpy as np
import pickle

import sys
sys.path.append('..') # make masonweaver_analytic available here

#%% Set parameters

SAVEFIG = False

#%% import matplotlib_werts_figsty1
from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['figure.figsize'] = (3.149, 3.5) # 80mm x 50mm
rcParams['figure.subplot.left'] = 0.21 
rcParams['figure.subplot.right'] = 0.94    
rcParams['figure.subplot.bottom'] = 0.12
rcParams['figure.subplot.top'] = 0.94

#%% Read results

file_name = 'FIG_IntS_ISS_comp_time_measure.pickle'

with open(file_name, 'rb') as f:
    y = pickle.load(f)
    tau = pickle.load(f)
    alpha = pickle.load(f)
    t_PWS = pickle.load(f)
    t_TISS = pickle.load(f)
    eps = pickle.load(f)
    tau_s = pickle.load(f)
    
#%% Visualization
plt.figure(1)
plt.clf()

plt.subplot(211)
plt.semilogy(tau, np.asarray(t_PWS), '-b', label=r'IntS/ISS switch')
plt.plot(tau, np.asarray(t_TISS),'-r', label=r'pure ISS')
plt.plot([tau_s, tau_s], [0, np.max(t_TISS),], ':g')
plt.legend(frameon=False)
plt.ylabel('calc. time [s]')

plt.subplot(212)
plt.semilogy(tau[0:],np.asarray(eps[0:]))
plt.plot([tau_s, tau_s], [0, np.max(eps[0:])*10], ':g')
plt.ylabel('MSE')
plt.xlabel('$\\tau$')

#plt.subplots_adjust(wspace=0.4, 
#                    hspace=0.4)

if SAVEFIG:
    plt.savefig('FIG_IntS_ISS_comp_time.png', dpi=450)
