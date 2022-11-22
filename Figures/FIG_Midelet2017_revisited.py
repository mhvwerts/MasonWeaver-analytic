# -*- coding: utf-8 -*-
"""
FIGURE. Take published experimental sedimentation profiles (Fig. 2 from 
Ref. [1]) and compared with theoretical profiles from the analytic Mason-Weaver
solution.


Reference
[1]  J. Midelet, A. H. El-Sagheer, T. Brown, A. G. Kanaras, M. H. V. Werts
     “The Sedimentation of Colloidal Nanoparticles in Solution and 
     Its Study Using Quantitative Digital Photography.”
     Part. Part. Syst. Charact. 2017, 34, 1700095.
     https://doi.org/10.1002/ppsc.201700095.
"""

SAVEFIG = False


#%% Imports

from pathlib import Path
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

# Import 'real-world' analytic Mason-Weaver profile function
import sys
sys.path.append('..') 
from masonweaver_analytic import MW_c_profile


#%% Configure

# Data files containing the experimental optical density (OD) vs height (z),
# traces obtained from digital photographs as explained in Ref. [1]
#

fdir = './Midelet2017_Fig2_exp_data/'
flist = ["DSCN0063_ODtrace_green.hdf5",
         "DSCN0064_ODtrace_green.hdf5",
         "DSCN0066_ODtrace_green.hdf5",
         "DSCN0068_ODtrace_green.hdf5",
         "DSCN0069_ODtrace_green.hdf5",
         "DSCN0070_ODtrace_green.hdf5",
         "DSCN0071_ODtrace_green.hdf5",
         "DSCN0072_ODtrace_green.hdf5",
         "DSCN0073_ODtrace_green.hdf5",
         "DSCN0074_ODtrace_green.hdf5",
         "DSCN0075_ODtrace_green.hdf5",
         "DSCN0077_ODtrace_green.hdf5",
         "DSCN0078_ODtrace_green.hdf5",
         "DSCN0079_ODtrace_green.hdf5",
         "DSCN0080_ODtrace_green.hdf5",
         "DSCN0081_ODtrace_green.hdf5"
         ]


#%% Process and calculate

# Read individual traces from files
tt = [] # time data (seconds since start of experiment)
traces = [] # individual OD vs z traces (dict with array)
for fname in flist:
    fp = Path(fdir, fname)
    with h5.File(fp, 'r') as f:
        print(fname, '\t', f.attrs['photo_timestamp_iso'])
        tt.append(f.attrs['seconds_from_exp_start'])
        z = np.array(f['z'])
        OD = np.array(f['OD_exp'])
        traces.append({'z_exp': z,
                       'OD_exp': OD})

# Time points
tt = np.array(tt)

# Create data matrix with all experimental traces 
z_exp = traces[0]['z_exp'] # all traces were sampled on the same z grid
OD_exp = np.zeros((len(z_exp), len(traces)))
for ix, trace in enumerate(traces):
    OD_exp[:, ix] = trace['OD_exp']

# Mason-Weaver fit parameters from original data (Ref. [1], Figure 2)
# Fit result (20160728)
p_opt = np.array([5.12583364e-12,   7.76175686e-09,   4.47386306e-01,
                  1.64675504e+01])
D = p_opt[0]
sg = p_opt[1]
OD0 = p_opt[2]
z_max_mm = p_opt[3]

# Create data matrix with all theoretical traces
Nz = 500
OD_theo = np.zeros((Nz, len(tt)))
for ix, t in enumerate(tt):
    z, c_ = MW_c_profile(t, z_max_mm/1000, D, sg, Nz=Nz)
    OD_theo[:, ix] = OD0 * c_
z_theo = z * 1000 # from meters to mm


#%% Create plot
# Plot settings for figures
from matplotlib import rcParams

rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['font.size'] = 9
rcParams['figure.figsize'] = (3.149, 2.8) # 80mm x 50mm
rcParams['figure.subplot.left'] = 0.15
rcParams['figure.subplot.right'] = 0.98    
rcParams['figure.subplot.bottom'] = 0.15
rcParams['figure.subplot.top'] = 0.98

# Create Figure
plt.figure(1)
plt.clf()

# Plot experimental data
plt.subplot(211)
plt.locator_params(axis='y',nbins=5)
plt.tick_params(direction='in', top=True, bottom=True, left=True, right=True)
plt.gca().axes.xaxis.set_ticklabels([])
plt.plot(z_exp, OD_exp,'r')
plt.arrow(12.5, 0.32, -5.9, -0.14,
          width = 0.010, head_width = 0.07, head_length = 0.6,
          fc='k', ec='k')
#plt.xlabel('z position / mm')
plt.ylabel('OD')
plt.text(8.3, 0.52, "exp. optical density")
plt.xlim(-0.8,18)
plt.ylim(-0.05,0.65)

# Plot analytic Mason-Weaver traces
plt.subplot(212)
plt.locator_params(axis='y',nbins=5)
plt.tick_params(direction='in', top=True, bottom=True, left=True, right=True)
plt.plot(z_theo, OD_theo, 'k')
plt.xlabel('z position / mm')
plt.ylabel('calc. OD')
plt.text(4.0, 0.52, "analytic Mason-Weaver soln.")
plt.xlim(-0.8,18)
plt.ylim(-0.05,0.65)

plt.savefig('FIG_Midelet2017_revisited.png',
            dpi=300)
plt.show()
