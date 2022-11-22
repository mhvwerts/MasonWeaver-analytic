#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FIGURE. Comparison between the curves plotted in Figure 1 of the original 
work by Mason & Weaver [1] and those obtained with the present algorithm.

[1] Max Mason, Warren Weaver,
    "The Settling of Small Particles in a Fluid",
    Phys. Rev. 1924, 23, 412
"""

DRAFT = True


#%% Import 
import csv
from matplotlib import pyplot as plt
import numpy as np

import sys
sys.path.append('..') # make masonweaver_analytic available here

from masonweaver_analytic import MW_adim


#%% Set formatting for Figure plots

from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['figure.figsize'] = (3.149, 2.362) # 80mm x 50mm
rcParams['figure.subplot.left'] = 0.19  
rcParams['figure.subplot.right'] = 0.94    
rcParams['figure.subplot.bottom'] = 0.19
rcParams['figure.subplot.top'] = 0.94


#%% Read and process the CSV file extracted from Fig. 1,  Mason&Weaver 1924
#   using the Engauge software

f = open('MasonWeaver1924_Fig1_digitized/Fig_1.csv', 'r')

raw_data = []
with f:
    reader = csv.reader(f, delimiter=",")
    for row in reader:
        raw_data.append(row)

x =[]
c0=[]
c1=[]
c2=[]
c3=[]
c4=[]
        
for i in range(1, len(raw_data)):    
    x.append(float(raw_data[i][0]))
    c0.append(float(raw_data[i][1]))  
    c1.append(float(raw_data[i][2]))
    c2.append(float(raw_data[i][3]))
    c3.append(float(raw_data[i][4]))
    c4.append(float(raw_data[i][5]))


#%% Set parameters to go with Mason&Weaver 1924 and present evaluation

ny = 200 # Define 1D domain with nx solution points
y = np.linspace(0.0, 1.0, ny)

alpha = 0.025 # Value for Figure 1 in M&W 1924


#%% Compute adimensional solution using MW_adim

c0_th = MW_adim(y, 0.05, alpha)
c1_th = MW_adim(y, 0.25, alpha)
c2_th = MW_adim(y, 0.5, alpha)
c3_th = MW_adim(y, 1., alpha)
c4_th = MW_adim(y, None, alpha) # tau = None evaluates stead-state solution


#%% Plot


plt.figure(1)
plt.clf()
ax = plt.gca()
plt.plot(x[:110],c0[:110],':b')
plt.plot(x,c1,':r')
plt.plot(x,c2,':y')
plt.plot(x,c3,':g')
plt.plot(x,c4,':k')
plt.plot(c0_th,y,'b')
plt.plot(c1_th,y,'r')
plt.plot(c2_th,y,'y')
plt.plot(c3_th,y,'g')
plt.plot(c4_th,y,'k')
plt.ylim(1,0)
plt.xlim(0,7)

# annotation
ax.text(2.0, 0.10, '$\\tau = 0.05$')
plt.plot([1.03, 1.97], [0.117, 0.08], 'k:', linewidth=0.8)
ax.text(2.4, 0.25, '$\\tau = 0.25$')
plt.plot([0.68, 2.37], [0.280, 0.235], 'k:', linewidth=0.8)
ax.text(2.7, 0.47, '$\\tau = 0.5$')
plt.plot([0.56, 2.63], [0.49, 0.45], 'k:', linewidth=0.8)
ax.text(2.9, 0.68, '$\\tau = 1.0$')
plt.plot([0.43, 2.85], [0.72, 0.66], 'k:', linewidth=0.8)
plt.plot([0.25, 2.85], [0.79, 0.66], 'k:', linewidth=0.8)
ax.text(3.0, 0.85, '$\\tau \\rightarrow \\infty$')
plt.plot([0.25, 2.97], [0.860, 0.83], 'k:', linewidth=0.8)

# plt.legend()
plt.xlabel(r'$c(y, \tau)/c_0$')
plt.ylabel(r'$y$')

if not DRAFT:
    plt.savefig('FIG_MW1924_Comparison.png', dpi=300)

plt.show()

