# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
from pyquantize.numerics.bases.fourier import basis as fourier_basis
from pyquantize.systematics.system_phase import System_Phase
from pyquantize.numerics.spaces.space_1D import Space_1D
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

# parameters for the qubit
cshunt_params={ 'E_j': 186e9, 'E_c': .98e9, 'alpha': .65, 'f': .498 }   
            
# define the flux qubit problem
space=Space_1D([-2*np.pi,2*np.pi],Nx=10000)
        
problem_func= lambda **params:\
    {'space': space,
    'mass': 1/(8*params['E_c']),
    'potential': -params['E_j']*(2*np.cos(space.x/2)+
        params['alpha']*np.cos(space.x-2*np.pi*params['f']))} 
        

cshunt=System_Phase(problem_func,fourier_basis,Nbasis=50,**cshunt_params)
cshunt.visualize_levels(5)
cshunt.visualize_state('lowest_four')
print np.diff(cshunt._sys_1d._energies[:5])