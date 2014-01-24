# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
from __future__ import division
from pyquantize.numerics.spaces.space_1D import Space_1D
from pyquantize.numerics.bases.fourier import basis as fourier_basis
from pyquantize.numerics.bases.coord_SHO import basis as sho_basis
from pyquantize.systematics.system_1d import System_1D
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

space=Space_1D([-3,3],'vanish',3000)
mass=2
def problem(freq):
    return {'space': space,
            'mass': mass,
            'potential': 1/2*mass*(2*np.pi*freq)**2*np.power(space.x,2) }
            
basis_func=fourier_basis

f=1
shof=System_1D(problem,basis_func,Nbasis=90,freq=f)
shof.visualize_levels(5)
shof.visualize_state('lowest_four')
shof.visualize_H(20)
print np.round(shof._energies[:5],4)
#
#basis_func=sho_basis
#
#sho=System_1D(problem,basis_func,freq=f)
#sho.visualize_levels(5)
#sho.visualize_state('lowest_four')
#sho.visualize_H(20)
#print np.round(sho._energies[:5],4)
