# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 

"""
from pyquantize.systematics.system import System
from pyquantize.systematics.system_1d import System_1D
import numpy as np
import matplotlib.pyplot as plt

class System_Phase(System):
    
    def __init__(self, problem_func, basis_func, **others):
        def modified_pf(**params):
            prob=problem_func(**params)
            prob.update(mass=prob['mass']/(4*np.pi**2))
            return prob
            
        self._sys_1d=System_1D(modified_pf,basis_func,**others)
        
    def set_params(self,**params):
        self._sys_1d.set_params(**params)
    
    def _eval_simple_op(self,operator):
        if operator=='H':
            return self._sys_1d._eval_simple_op('H')
        if operator=='delta':
            return self._sys_1d._eval_simple_op('x')
        elif operator=='n':
            return self._sys_1d._eval_simple_op('p')*2*np.pi
        else:
            return System._eval_simple_op(self,operator)         
            
    def visualize_levels(self,n):
        self._sys_1d.visualize_levels(n)
        plt.xlabel('$\delta$')

    def visualize_state(self,n):
        self._sys_1d.visualize_state(n)
        
    def _plot_state(self,n,ax1):
        self._sys_1d._plot_state(n,ax1)
        plt.xlabel('$\delta$')

    def visualize_H(self,n):
        self._sys_1d.visualize_H(n)
        
    def get_energies(self):
        return self._sys_1d.get_energies()
        
    def get_dimension(self):
        return self._sys_1d.get_dimension()