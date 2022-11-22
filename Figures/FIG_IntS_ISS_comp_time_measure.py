#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Measure the computation time needed to evaluate the IntS and ISS to the
dimensionless Mason-Weaver equation
"""
#%%
import timeit
import numpy as np
import pickle


#%% Import MWE analytical lib.
import sys
sys.path.append('..') # make masonweaver_analytic available here

from masonweaver_analytic import MW_adim_ISS, MW_adim
from masonweaver_analytic import tau_switch, MSE



#%% Calculation parameters (some need to be specified twice)

mysetup = '''
import numpy as np
from masonweaver_analytic import MW_adim_ISS, MW_adim

#%% Adimensional space vector
ny = 200 # Define 1D domain with ny solution points
y = np.linspace(0.0,1.0,ny) # Generate adimensionnal vector

#%% Péclet number
alpha = 0.02 # Inverse of the Péclet number [.]
'''

# Adimensional time vector
steps = 250 # Define number of steps [.]
tau = np.linspace(0.001,1,steps+1)

# Adimensional height vector
ny = 200 # Define 1D domain with ny solution points
y = np.linspace(0.0,1.0,ny) # Generate adimensional vector

# Alpha 
alpha = 0.02 # Inverse of the Péclet number [.]


#%% Corresponding tau_switch
tau_s =  tau_switch(alpha)
print('Switching at tau = ', tau_s)



#%% Timeit statement 
Nruns = 100

t_PWS = []
for i in range(0, tau.size):
    print(i,'/', tau.size)
    mycode = '''c_ = MW_adim(y,'''+str(tau[i])+''', alpha)'''
    t_PWS.append(timeit.timeit(setup = mysetup,
                     stmt = mycode, number = Nruns)/Nruns)
    
t_ISS = []
for i in range(0, tau.size):
    print(i,'/', tau.size)
    mycode = '''c_ = MW_adim_ISS(y, '''+str(tau[i])+''', alpha)'''
    t_ISS.append(timeit.timeit(setup = mysetup,
                     stmt = mycode, number = Nruns)/Nruns)
    
#%% Computation of the MSE error
eps = []
for i in range(0,tau.size):
    eps.append(MSE(MW_adim(y, tau[i], alpha),MW_adim_ISS(y, tau[i], alpha)))
    
#%% Save concentration profile as pickle

file_name = 'FIG_IntS_ISS_comp_time_measure.pickle'
with open(file_name, 'wb') as f:
    # Pickle the 'data' dictionary 
    pickle.dump(y, f)
    pickle.dump(tau, f)
    pickle.dump(alpha, f)
    pickle.dump(t_PWS,f)
    pickle.dump(t_ISS,f)
    pickle.dump(eps,f)
    pickle.dump(tau_s,f)
