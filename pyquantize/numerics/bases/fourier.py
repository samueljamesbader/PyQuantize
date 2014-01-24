# -*- coding: utf-8 -*-
"""
Module: qsim.numerics.bases.fourier
Author: Sam Bader

Sinusoidal basis functions.

"""
from __future__ import division
import numpy as np

def basis(problem, Nbasis, tolerance=.1):
    ''' Returns a basis of sinusoidal functions.  Each column of the matrix is
    one basis function (one k-value) evaluated over the x values of the space.
    
    space: the space for the wavefunctions
    Nbasis: the number of basis functions to return
    potential: 1-D array-like potential [not used for this basis function]
    boundary: what boundary conditions to obey: either 'vanish' for sinusoids
    which vanish at the boundary, or 'PBC' for periodic boundary conditions.
    
    '''
    space=problem['space']
    
    if space.boundary == 'vanish':
        # Condition to vanish on the boundary
        n=np.matrix(np.arange(0,Nbasis))
        k=(np.pi/space.L)*(n+1).T
        
        # Smallest length scale is
        lambda_min = 2*np.pi/k[-1,-1]        
        
        # Normalized sin basis centered at left boundary
        basis=np.sqrt(2/space.Nx)*np.sin((space.x-space.xmin)*k.T)

    elif space.boundary == 'PBC':
        
        # Half the states [or less if Nbasis is odd] should be sines
        n_sin=np.matrix(np.arange(0,np.floor(Nbasis/2)))
        k_sin=(2*np.pi/space.L)*(n_sin+1).T
        
        # The other half [or more if Nbasis is odd] should be cosines
        # Note that a k=0 is included for cosines, which is why cosines should
        # naturally take the "one more" if Nbasis is odd.
        n_cos=np.matrix(np.arange(0,np.ceil(Nbasis/2)))
        k_cos=(2*np.pi/space.L)*(n_cos).T
        
        # Normalized wavefunctions.  Note the patch repairing the normalization
        # of the k=0 state.
        sins=np.sqrt(2/space.Nx)*np.sin((space.x-space.xmin)*k_sin.T)
        coss=np.sqrt(2/space.Nx)*np.cos((space.x-space.xmin)*k_cos.T)
        coss[:,0]=coss[:,0]/np.sqrt(2);
        
         # Smallest length scale is
        lambda_min = 2*np.pi/k_cos[-1,-1]        
        
        # Compile them altogether
        basis=np.matrix(np.zeros([space.Nx,Nbasis]))
        basis[:,::2]=coss;
        basis[:,1::2]=sins;
    
    max_dx=lambda_min*tolerance
    if max_dx < space.dx:
        print "BASIS CANNOT BE FULLY RESOLVED WITH CURRENT dx"
        print\
        "Maximum dx for this basis with this tolerance is {}".format(max_dx)
        
    return basis