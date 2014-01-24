# -*- coding: utf-8 -*-
"""
Module qsim.numerics.bases.SHO
Author: Sam Bader

Simple Harmonic Oscillator basis functions in a position space.

"""

from __future__ import division
import numpy as np
import qsim.numerics.util as util
import numpy.polynomial.hermite as herm
import scipy.linalg as lin


def H(x,n):
    '''Returns a matrix of Hermite polynomials.  Each column represents one
    H_n(x) polynomial (physicist's convention) evaluated over a range of x's.
    
    x: 1-D array-like of x-values at which to evaluate the polynomials
    n: 1-D array-like of integers indicating which polynomial to evaluate
    
    '''

    # Convert the listing of polynomial numbers into a matrix of columns where
    # each column has a 1 marking the polynomial to use and zeros elsewhere
    # (because this is what the numpy hermval asks for).
    nvals=util.force_array(n)
    c=np.zeros([np.max(nvals)+1,len(nvals)])
    for i,n_i in enumerate(nvals):
        c[n_i,i]=1
    
    # Evaluate
    return herm.hermval(util.force_array(x),c).T


def phi_n(space,x0,n,l):
    ''' Returns a matrix of unnormalized SHO wavefunctions centered about x0.
    
    space: the Space for the wavefunctions
    x0: the center of the potential well
    n: 1-D array-like of integers indicating which level (indexed from 0)
    l: length scale sqrt(1/m*omega) in the Hamiltonian (note h=1)
    
    The wavefunctions are not normalized here, because the standard
    normalization for SHO wavefunctions would fail (because these
    wavefunctions are given on only a finite domain).  The function which uses
    these wavefunctions to form a basis will thus need to apply Gram-Schmidt.
    
    '''
    
    # Center at x0
    x=space.x-x0   
    
    # n-dependent prefactor for normalization (single-row matrix)
#    n=np.matrix(n)
#    n_prefactor=\
#        (np.sqrt(space.L/(space.Nx*np.multiply(
#            np.power(2,n),
#            sp.misc.factorial(n))))\
#            *q**(1/8)/np.pi**(1/4))
    
    # Exponential suppression for large x (single-column matrix)
    exp_supp=np.exp(-np.power(x/l,2)/2)
    
    # Hermite polynomial (Nx-by-size(n) matrix)
    herm_mod=H(x/l,n)
    
    # Multiply elementwise
#    return np.multiply(n_prefactor,np.multiply(exp_supp,herm_mod))
    return np.multiply(exp_supp,herm_mod)


def basis(problem,Nbasis):
    ''' Returns a basis of SHO wavefunctions.  Each column of the matrix is one
    basis function evaluated over the x values of the space.
    
    The potential is examined for a minimum, and the wavefunctions are the SHO
    basis states about that minima, but re-orthonormalized in this space.
    
    space: the space for the wavefunctions
    Nbasis: the number of basis functions to return
    potential: 1-D array-like potential
    
    '''
    
    space=problem['space']
    potential=problem['potential']
    mass=problem['mass']
    
    # Find the minima and examine the curvature there to extract the q-value
    # as defined in the documentation of phi_n()
    n=np.arange(Nbasis)
    min_index=np.argmin(potential)
    x0=space.x[min_index,0]
    l=1/(4*np.pi**2*mass*(space.D*space.D*potential)[min_index,0])**.25
    print "l={}".format(l)
    
    # Get the basis states, apply Gram-Schmidt
    exact_SHO=phi_n(space,x0,n,l)
    re_orthogonalized,_ = lin.qr(exact_SHO,mode='economic')
    
    print "Automatic Checking of dx not yet implemented for SHO basis"
    
    return np.matrix(re_orthogonalized)