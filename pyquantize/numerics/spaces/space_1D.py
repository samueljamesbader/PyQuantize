# -*- coding: utf-8 -*-
"""
Module: numpy.spaces.Space
Author: Sam Bader

Define a 1-D discretized space
"""
from __future__ import division
import numpy as np
import scipy.sparse as spar
from pyquantize.numerics.util import force_array

class Space_1D:
    ''' Collects all the properties for a space into one object.
    
    Recorded properties [for convenience]
    --
    xmin: lower boundary of domain
    xmax: upper boundary of domain
    L: length of domain
    Nx: number of discretized chunks
    dx: width of discretized chunk
    
    Useful variables
    --
    x: single column matrix of x-values for the space
    D: sparse matrix representing the derivative operator.
    
    '''
    
    def __init__(self,domain,boundary='PBC',Nx=1000):
        ''' Define a 1-D space over a domain.
        
        domain: tuple of (xmin, xmax)
        Nx: integer number of chunks to discretize space
        
        '''

        # Record properties of the space
        self.xmin=domain[0]
        self.xmax=domain[1]
        self.L=self.xmax-self.xmin
        self.Nx=Nx
        self.boundary=boundary
        
        # Discretize
        self.x=np.matrix(np.linspace(self.xmin,self.xmax,self.Nx,endpoint=False)).T
        self.dx=self.L/self.Nx
        
        # Define derivative matrix
        self.D=spar.diags([np.ones(Nx-1),-np.ones(Nx-1),np.array([-1]),
            np.array([1])],[1,-1,Nx-1,-(Nx-1)])/(2*self.dx)
            
        self.q=spar.diags([force_array(self.x)],[0])
        self.p=-1j*self.D/(2*np.pi)
