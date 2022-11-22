#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
masonweaver_analytic.py :

Numerical evaluation of analytic solutions of the Mason-Weaver equation

    Lancelot Barthe, 
    École normale supérieure de Rennes
    CNRS, lab. SATIE
    Campus de Ker Lann, F-35170 Bruz, France

    Martinus H. V. Werts, 
    École normale supérieure de Rennes 
    CNRS, lab. SATIE,
    Campus de Ker Lann, F-35170 Bruz, France


Distributed under the GPL-compatible CeCILL license.
Read the license text at the end of this file before using this software.


References :

[1] Max Mason, Warren Weaver,
    "The Settling of Small Particles in a Fluid",
    Phys. Rev. 1924, 23, 412

'''

#%% Imports

import numpy as np
from scipy.special import erf



#%% Globals

# Definition of ranges of validity of different solutions.
# See accompanying paper to see how these were established.
#
# These parameters were optimized for Python-numpy-scipy on Intel processor.
# They are expected to be suitable also for MATLAB, and other standard
# systems for scientific floating point calculations (FORTRAN, C,
# JavaScript etc.)
#
ALPHA_SWITCH1 = 0.0015  # frontier between hiPe and hiPe2
ALPHA_SWITCH2 = 0.02    # begin IntS/ISS switch zone
ALPHA_SWITCH3 = 20.0    # begin pure ISS zone

TAU_SWITCH_A1 = 8.134    # parameters for calculation of tau_switch
TAU_SWITCH_K1 = 0.9666
TAU_SWITCH_A2 = 0.02934
TAU_SWITCH_K2 = -0.8634
TAU_SWITCH_ALPHA_MAX = 0.04625 # exp((ln(A2)-ln(A1))/(k1-k2))


# Constants
MAX_ISS_TERMS = 100000  # Maximum number of ISS terms to be evaluated
                        # in MW_adim_ISS(). Raise error if number of terms 
                        # exceed this maximum.
                        
# Diagnostics
ISS_Nterms_latest = -1  # Number of terms contributing to the latest evaluation
                        # of MW_adim_ISS()
                        


#%% Function definitions

def MW_adim_steadystate(y, alpha):
    """
    Steady-state sedimentation profile
    

    Parameters
    ----------
    y : numpy vector or float
        Adimensional spatial coordinates (depth).
    alpha : float
        Inverse of the Péclet number.

    Returns
    -------
    c_ss : numpy vector
        Concentration profile at steady state.
    """
    # Steady-state solution
    c_num = np.exp(y/alpha) # numerator
    c_denom = (alpha * (np.exp(1/alpha)-1)) # denominator
    c_ss = c_num / c_denom # full expression
    return c_ss
    


def MW_adim_ISS(y, tau, alpha):
    '''
    Compute the solution of the dimensionless MWE using the Infinite Series
    Solution (ISS)

    Due to floating-point limitations, this routine fails to give reliable 
    results when alpha < 0.02. (64 bits floating point).
    
    For a certain range of (alpha, tau), the alternative MW_adim_IntS will
    calculate a suitable solution much faster. See accompanying paper.
    
    
    Parameters
    ----------
    y : numpy vector or float
        Adimensional spatial vector [.].
    tau : float
          Adimensional time [.].
    alpha : float
            Inverse of the Péclet number [.].

    Returns
    -------
    c_ : numpy vector of length len(y)
        Concentration profile along y-axis calculated
        for a given alpha and tau.

    Continues adding terms until the latest term numerically evaluates as zero.
    (See Main Text). The number of terms is stored in the global variable
    `ISS_Nterms_latest`

    As a precaution, the maximum number of terms to be evaluated is set by the
    global variable `MAX_ISS_TERMS`. An error is raised if the number of terms
    exceeds this value. If, by luck, this occurs, please contact M.H.V. Werts!
    '''
    global MAX_ISS_TERMS
    global ISS_Nterms_latest

    # Steady-state solution
    c_ss = MW_adim_steadystate(y, alpha)

    # Transient solution (series)
    ## Prefactor of the series
    s_factor = 16 * alpha**2 * np.pi * np.exp((2*y - tau)/(4*alpha))
    ## Initialise the sum over series terms (s_m)
    S_m = np.zeros_like(y) # 
    for m in range(1, MAX_ISS_TERMS+1):
        # Decomposition of the series terms in four subterms 
        # to simplify the expressions
        s_a_m = np.exp(-alpha * m**2 * np.pi**2 * tau)
        s_b_m = m * (1 - (-1)**m * np.exp(-1/(2*alpha)))
        s_c_m = (np.sin(m*np.pi*y) +\
                 2.0*np.pi*m*alpha*np.cos(m*np.pi*y))
        s_d_m = (1 + 4 * np.pi**2 * m**2 * alpha**2)**2
        if s_a_m/s_d_m == 0: # term vanishes numerically, end loop
            ISS_Nterms_latest = m-1 # store number of terms used
            break
        s_m = ((s_a_m*s_b_m*s_c_m)/s_d_m)
        S_m += s_m
    if m >= MAX_ISS_TERMS:
        raise Exception('Maximum number of ISS terms'\
                        ' exceeded in MW_adim_ISS().')
    c_t = s_factor*S_m # Compute transient solution
    c_ = c_ss + c_t # Compute overall solution
    return c_



def MW_adim_IntS(y, tau, alpha):
    '''
    Compute the integral solution (IntS) of the dimensionless.

    Faster computation of the solution compared to ISS for a certain range
    (alpha, tau). Limited domain of validity, see accompanying paper.


    Parameters
    ----------
    y : numpy vector
        Adimensional vector spatial vector [.].
    tau : float
          Adimensional time [.].
    alpha : float
            Inverse of the Péclet number [.].

    Returns
    -------
    c_ : numpy vector of length len(y)
        Concentration profile along y-axis calculated
        for a given alpha and tau.

    '''
    # First term of the expression
    a_1 = (1/2)*np.exp((y-1)/alpha)
    a_2 = 1 + ((y+tau-1)/alpha)
    a_31 = erf((2-y-tau)/np.sqrt(4*alpha*tau))
    a_32 = erf((1-y-tau)/np.sqrt(4*alpha*tau))
    A   = a_1*a_2*(a_31-a_32)

    # Second term
    b_1  = np.sqrt(tau/(alpha*np.pi))*np.exp((y-1)/alpha)
    b_21 = np.exp(-((2-y-tau)/np.sqrt(4*alpha*tau))**2)
    b_22 = np.exp(-((1-y-tau)/np.sqrt(4*alpha*tau))**2)
    B   = b_1*(b_21 - b_22)

    # Third term
    c_1 = (1/2)*np.exp(y/alpha,dtype = np.float64)
    c_2 = 1 + ((y+tau)/alpha)
    c_31 = erf((1+y+tau)/np.sqrt(4*alpha*tau))
    c_32 = erf((y+tau)/np.sqrt(4*alpha*tau))
    C   = c_1*c_2*(c_31-c_32)

    # Fourth term
    d_1 = np.sqrt(tau/(alpha*np.pi), dtype=np.float64)\
          * np.exp(y/alpha, dtype = np.float64)
    d_21 = np.exp(-((1+y+tau)\
           / np.sqrt(4*alpha*tau, dtype=np.float64))**2, dtype=np.float64)
    d_22 = np.exp(-((y+tau)\
           / np.sqrt(4*alpha*tau, dtype=np.float64))**2, dtype = np.float64)
    D   = d_1*(d_21-d_22)

    # Fifth term
    e_11 = erf((1-y+tau)/np.sqrt(4*alpha*tau))
    e_12 = erf((tau-y)/np.sqrt(4*alpha*tau))
    E = (1/2)*( e_11 - e_12 )

    # Combination of all terms
    c_ = A - B + C + D + E

    return c_



def MW_adim_hiPe(y, tau, alpha):
    """Evaluate analytic expression describing Mason-Weaver sedimentation
    profiles at higher Péclet numbers (low values for alpha) - "hiPe"
    
    Valid for 0.0015 < alpha < 0.02.
    

    Parameters
    ----------
    y : numpy vector (1D ndarray)
        Vector containing adimensional height coordinates.
    tau : float
        Adimensional time.
    alpha : float
        Adimensional parameter, inverse of Péclet number.

    Returns
    -------
    c_ : numpy vector (1D ndarray)
        Vector containing relative concentration as a function of `y`, 
        calculated at time `tau` for a system defined by `alpha`.

    """
    A = 1 + erf((y-tau)/np.sqrt(4*alpha*tau))

    KA1 = np.sqrt((4*alpha*tau)/np.pi)
    KA2 = np.exp(-(tau-1)**2/(4*alpha*tau)) - np.exp(-tau/(4*alpha))
    KB = (tau-1)*erf((1-tau)/np.sqrt(4*alpha*tau))
    KC = tau*erf(tau/np.sqrt(4*alpha*tau))
    K = 1 - KA1*KA2 + KB + KC

    B = K * np.exp(y/alpha)/(alpha*(np.exp(1/alpha)-1))

    c_ = 0.5*(A + B)
    return c_



def MW_adim_hiPe2(y, tau, alpha):
    """Evaluate analytic expression describing Mason-Weaver sedimentation
    profiles at the highest Péclet numbers (low values for alpha) - "hiPe2"
        
    Valid for 0 < alpha < 0.0015. 
    

    Parameters
    ----------
    y : numpy vector (1D ndarray)
        Vector containing adimensional height coordinates.
    tau : float
        Adimensional time.
    alpha : float
        Adimensional parameter, inverse of Péclet number.

    Raises
    ------
    Exception
        This function is only valid for 0 <= y <= 0.99. Raises an Exception
        if this vector contains elements that are out of bounds.

    Returns
    -------
    c_ : numpy vector (1D ndarray)
        Vector containing relative concentration as a function of `y`, 
        calculated at time `tau` for a system defined by `alpha`.

    """
    if not ((min(y)>=0) and (max(y)<=0.99)):
        raise Exception('y out of domain of validity of MW_adim_hiPe2')

    A = 1 + erf((y-tau)/np.sqrt(4*alpha*tau))

    c_ = 0.5*A
    
    return c_




def MW_adim(y, tau, alpha):
    '''
    Compute the piece-wise analytic solution of the dimensionless MWE
    

    Parameters
    ----------
    y : numpy vector
        Adimensional vector spatial vector [.].
    tau : float or None
          Adimensional time [.].
          If tau == None, return the steady state profile
    alpha : float
            Inverse of the Péclet number [.].

    Returns
    -------
    c_ : Numpy Vector of size ny
        Concentration profile along y-axis estimated
        for a given alpha and tau.

    '''
    # Check if parameters within domain of validity

    if alpha <= 0:
        raise Exception('alpha should be positive.')
    if (min(y) < 0.0) or (max(y) > 1.0):
        raise Exception('y out of range')
    # Choose appropriate analytic solution to be evaluated
    if tau is None:
        c_ = MW_adim_steadystate(y, alpha)
    else:
        if tau < 0:
            raise Exception('tau should be non-negative.')
        elif tau == 0:
            c_ = np.ones_like(y)
        else:
            if alpha < ALPHA_SWITCH1:
                c_ = MW_adim_hiPe2(y, tau, alpha)
            elif ALPHA_SWITCH1 <= alpha < ALPHA_SWITCH2:
                c_ = MW_adim_hiPe(y, tau, alpha)
            elif ALPHA_SWITCH2 <= alpha < ALPHA_SWITCH3:
                tau_sw = tau_switch(alpha)
                if tau < tau_sw:
                    c_ = MW_adim_IntS(y, tau, alpha)      
                else:
                    c_ = MW_adim_ISS(y, tau, alpha)
            else:
                c_ = MW_adim_ISS(y, tau, alpha)
    return c_



def tau_switch(alpha) :
    '''Compute tau_switch at which the PWS will switch between the
    IntS and the ISS

    Parameters
    ----------
    alpha : float
        Inverse of the Péclet number [.].
            

    Returns
    -------
    float
        Dimensionless tau for switching between ISS and IntS.

    '''
    if alpha < TAU_SWITCH_ALPHA_MAX :
        tau_switch = TAU_SWITCH_A1 * pow(alpha, TAU_SWITCH_K1)
    else : 
        tau_switch = TAU_SWITCH_A2 * pow(alpha, TAU_SWITCH_K2)

    return tau_switch 



def MW_c_profile(t, z_max, D, sg, Nz = 200, full_output = False):
    """Calculate MW concentration profile, 'real-world' (non-adimensional) 
    

    Parameters
    ----------
    t : float [seconds] or None
        Time for which the profile should be evaluated. If None, evaluate the
        steady-state (sedimentation equilibrium) solution (lim t -> infinity)
    z_max : float [meters]
        Height of the liquid column.
    D : float [m^2 s^-1]
        Diffusion coefficient of the nanoparticles.
    sg : float [m s^-1]
        Product of the sedimentation coefficient and the (gravitational
        acceleration.
    Nz : int, optional
        Number of (height) sampling the concentration profile.
        The default is 200.
    full_output : bool, optional
        If True, additionally return a dict containing all parameters,
        adimensional, and dimensional, describing the MW system.
        The default is False.

    Returns
    -------
    z : 1D numpy.ndarray of floats 
        Vector containing the height coordinates (bottom is at z = 0).
    c_ : 1D numpy.ndarray of floats
        Concentration of particles at (z, t), relative to concentration at
        the start of sedimentation (c/c_0)
    params: dict
        Parameters, adimensional, and dimensional, describing the MW system and
        used for the evaluation. Only return if `full_output == True`
    """
    z_0 = D/sg
    alpha = z_0/z_max
    t_sed = z_max/sg
    if t is None: # use None to calculate steady-state/equilibrium solution
        tau = None
    else:
        if t < 0:
            raise Exception('t should be non-negative.')
        tau = t/t_sed
    if alpha < ALPHA_SWITCH1:
        y = np.linspace(0, 0.99, Nz) # hiPe2
    else:
        y = np.linspace(0, 1.0, Nz)
    z = (1-y)*z_max
    c_ = MW_adim(y, tau, alpha)
    if full_output:
        params = {'y': y,
                  'tau': tau,
                  'alpha': alpha,
                  'z_0': z_0,
                  't_sed': t_sed,
                  't': t,
                  'z_max': z_max,
                  'D': D,
                  'sg': sg,
                  'Nz': Nz}
        return z, c_, params
    else:
        return z, c_



#%% Utility functions


def MSE(c_opt_, c_):
    """Compute the mean-squared error (MSE) between two concentration vectors

    Parameters
    ----------
    c_opt_ : numpy vector
        Tested concentration profile as a function of space.
    c_ : numpy vector
        Reference concentration profile (ground truth) as a function of space.

    Returns
    -------
    eps : float
        Mean squared error corresponding to the two concentration profiles.

    """
    eps = np.mean((c_opt_ - c_)**2)
    return eps



def ismonotone(cc, rtol=0.99, atol=1e-7):
    '''
    Determine if a given vector (set of discrete points containing function
    values) is monotonic (non-decreasing) within certain tolerances tolerances

    Parameters
    ----------
    cc : list
        Set of .
    rtol : float, optional
        Relative tolerance in [0,1]. The default is 0.99.
    atol : float, optional
        Absolute tolerance. The default is 1e-7. This 'absolute tolerance'
        exists to avoid problems when the values are very small


    Returns
    -------
    monotone : bool
        Boolean flag representing the monotony.
    i_err : int
        Index where the monotony is not respected anymore.

    '''
    i_err =-1
    monotone = True

    c_ = cc[0]
    for i, c in enumerate(cc[1:]):
        if c < 0:
            raise Exception('ismonotone only works with positive values')
        if c > c_:
            c_ = c
        elif (c+atol) < (rtol*c_):
            monotone = False
            i_err = i+1
            break
        else:
            pass

    if monotone:
        i_err = -1

    return monotone, i_err



#
# Copyright 2022, L. Barthe (ENS Rennes) and M.H.V. Werts (CNRS)

# lancelot.barthe at ens-rennes dot fr
# martinus.werts at ens-rennes dot fr

# This software is a computer program whose purpose is to evaluate numerically
# analytic solutions of the Mason-Weaver equation which describes sedimentation
# of small particles in a fluid, taking into account their Brownian diffusion.

# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". See also the "LICENSE" file supplied with this
# software package.

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
