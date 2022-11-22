# -*- coding: utf-8 -*-
"""
Retrieve the (alpha, tau) boundaries for each of the analytic solutions 
forming the PWS of the Mason-Weaver equation.
"""
SAVEFIG = False

#%% Imports

import numpy as np
from numpy import log10
import matplotlib.pyplot as plt

import sys
sys.path.append('..') # make masonweaver_analytic available here

# retrieve boundaries for different analytic solutions
from masonweaver_analytic import tau_switch
from masonweaver_analytic import ALPHA_SWITCH1, ALPHA_SWITCH2, ALPHA_SWITCH3

#%% Configure plotting

from matplotlib import rcParams
rcParams['lines.linewidth'] = 1.0
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['axes.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['figure.figsize'] = (3.149, 3.149*0.7) # 80mm x xxx mm
rcParams['figure.subplot.left'] = 0.19  
rcParams['figure.subplot.right'] = 0.9   
rcParams['figure.subplot.bottom'] = 0.19
rcParams['figure.subplot.top'] = 0.94
rcParams['legend.frameon'] = False

#%% Calculate boundaries

tautop = 1e2
taubot = 1e-4

hiPe_left_alpha = [ALPHA_SWITCH1, ALPHA_SWITCH1]
hiPe_left_tau = [taubot, tautop]

hiPe_right_alpha = [ALPHA_SWITCH2, ALPHA_SWITCH2]
hiPe_right_tau = [taubot, tautop]

tswitch_alphas = np.logspace(log10(ALPHA_SWITCH2), log10(ALPHA_SWITCH3), 1000)
tswitch_boundary = [tau_switch(alpha) for alpha in tswitch_alphas]

pureISS_alpha = [ALPHA_SWITCH3, ALPHA_SWITCH3]
pureISS_tau = [taubot, tswitch_boundary[-1]]




#%% Plot and annotate

style = 'k'

plt.figure(1)
plt.clf()
plt.loglog(hiPe_left_alpha, hiPe_left_tau, style)
plt.loglog(hiPe_right_alpha, hiPe_right_tau, style)
plt.loglog(tswitch_alphas, tswitch_boundary, style)
plt.loglog(pureISS_alpha, pureISS_tau, style)

plt.text(2.5e-4, 5e-2, 'hiPe2', rotation = 90)
plt.text(3.5e-3, 7e-2, 'hiPe', rotation = 90)
plt.text(1e-1, 1e-3, 'IntS')
plt.text(1e1, 1e-1, 'ISS')

plt.xlim(1e-4, 1e3)
plt.ylim(1e-4, 1e2)
plt.xlabel('$\\alpha$')
plt.ylabel('$\\tau$')

if SAVEFIG:
    plt.savefig('FIG_PWS_ISS_IntS_hiPe_hiPe2_boundaries.png',
                dpi = 300)
    
plt.show()
