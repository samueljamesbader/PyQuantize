# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""

from __future__ import division
import scipy.sparse as spar
import scipy.linalg as lin
import numpy as np
from pyquantize.numerics.util import force_array


def solve_tise(problem, basis, Nimportant, tolerance):
        kin_term=-problem['space'].D*problem['space'].D\
            /(8*np.pi**2*problem['mass'])
        pot_term=\
            spar.diags([force_array(problem['potential'])],[0])
            
        H=kin_term+pot_term        
        H_mat=basis.T*H*basis
        energies,eigvecs=lin.eigh(H_mat)
        eigvecs_in_x=basis*eigvecs
        
        inc_err=np.max(np.sqrt(np.sum(np.power(\
            np.abs(np.divide(H*eigvecs_in_x[:,:Nimportant],\
            energies[:Nimportant])\
            -eigvecs_in_x[:,:Nimportant]),2),0)))
            
        if(inc_err>tolerance):
            print "INCOMPLETENESS ERROR EXCEEDS TOLERANCE"
            
        
        rel_basis=(n for n in reversed(range(basis.shape[1]))\
            if np.max(eigvecs[n,:Nimportant])>tolerance).next()
        
        return energies,eigvecs_in_x,H_mat,inc_err