# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:14:49 2014

@author: sam
"""
from __future__ import division
from pyquantize.systematics.system import System
#import matplotlib.pyplot as plt
import numpy as np

class System_SHO(System):
    
    def __init__(self, Nbasis=50, **params):
            
        System.__init__(self,**params)
        self._Nbasis=Nbasis
        self._params_updated()
        
    def _params_updated(self):
        System._params_updated(self)
        self._energies= self._params['freq']*(np.arange(self._Nbasis)+1/2)        
    
    def _eval_simple_op(self,operator):
        
        N=self.get_dimension()        
        if operator=='H':
            return np.matrix(np.diag(self._energies[:N]))
        if operator=='a':
            return np.matrix(np.roll(np.diag(np.sqrt(np.arange(N))),-1,0))
        if operator=='a^T':
            return np.matrix(np.roll(np.diag(np.sqrt(np.arange(N))),-1,1))
        else:
            return System._eval_simple_op(self,operator)
    
    def get_dimension(self):
        return self._Nbasis
    
    def change_basis(self,basis_func):
        pass