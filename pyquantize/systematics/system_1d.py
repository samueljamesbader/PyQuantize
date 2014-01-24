# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:14:49 2014

@author: sam
"""
from pyquantize.systematics.system import System
from pyquantize.numerics.solvers.tise import solve_tise
import matplotlib.pyplot as plt
import numpy as np

class System_1D(System):
    
    def __init__(self, problem_func, basis_func,\
        Nimportant=5, tolerance=.01, Nbasis=50, **params):
            
        System.__init__(self,**params)
        self._problem_func=problem_func
        self._basis_func=basis_func
        self._Nimportant=Nimportant
        self._tolerance=tolerance
        self._Nbasis=Nbasis
        self._params_updated()
        
    def _params_updated(self):
        System._params_updated(self)
        self._problem=self._problem_func(**self._params)
        basis=self._basis_func(self._problem,self._Nbasis,self._tolerance)
        self._energies,self._eigenvecs,self._H_mat,self._inc_err=\
            solve_tise(self._problem,basis,self._Nimportant,self._tolerance)        
    
    def _eval_simple_op(self,operator):
        N=self._Nimportant
        if operator=='H':
            return np.diag(self._energies[:N])
        if operator=='x':
            return \
                self._eigenvecs[:,:N].T*\
                self._problem['space'].q*\
                self._eigenvecs[:,:N]
        elif operator=='p':
            return -.5j*\
                self._eigenvecs[:,:N].T*\
                self._problem['space'].D*\
                self._eigenvecs[:,:N]
        else:
            return System._eval_simple_op(self,operator)
    
    def visualize_levels(self,n):
        
        space=self._problem['space']
        potential=self._problem['potential']
        plt.figure()
        plt.plot(space.x,potential)
        for l in range(n):
            e=self._energies[l]
            level=[(e if potential[i,0]<e else None)\
                for i in range(space.Nx)]
            plt.plot(space.x,level,'r')
        plt.axis('tight')
        plt.xlabel('x')
        plt.ylabel('Energy')
        plt.title('Energy Levels')
    

    def visualize_state(self,n):
        plt.figure()
        
        if n=='lowest_four':
            for i in range(4):
                ax1=plt.subplot(2,2,i+1)
                self._plot_state(i,ax1)
            plt.tight_layout()
        else:
            ax1=plt.subplot(111);        
            self._plot_state(n,ax1)
        
    def _plot_state(self,n,ax1):


        space=self._problem['space']
        potential=self._problem['potential']        
        # Plot the potential, store the axis limits
        plt.plot(space.x,potential,
                 label='potential')
        plt.axis('tight')

        # Plot the energy level [only where it is above the potential]
        e=self._energies[n]
        level=[(e if potential[i,0]<e else None)\
            for i in range(space.Nx)]
        plt.plot(space.x,level,'r',
                 label='energy level {}'.format(n))
                 
        # Label the plot so far
        plt.xlabel('x')
        plt.ylabel('Energy')
        plt.title('State {}'.format(n))

        # Create a wavefunction overlay
        ax2=ax1.twinx()
        
        # Plot the wavefunction
        plt.plot(space.x,self._eigenvecs[:,n],'g',
                 label='state {}'.format(n))
        plt.ylabel('Wavefunction [au]')
        
        # Match the colors
        ax1.spines['right'].set_color('green')
        ax1.spines['top'].set_color('green')
        ax1.spines['left'].set_color('blue')
        ax1.spines['bottom'].set_color('blue')
        ax1.tick_params(axis='y', colors='blue')
        ax2.tick_params(axis='y', colors='green')        
        ax1.yaxis.label.set_color('blue')
        ax2.yaxis.label.set_color('green')
        
        # Re-tighten x-axis
        plt.axis('tight')

    def visualize_H(self,n):        
        plt.figure()
        plt.pcolormesh(np.flipud(np.asarray(np.abs(self._H_mat[:n,:n]))))
        plt.colorbar()
        plt.show()
    
    def get_energies(self,n):
        return self._energies[:n]
        
    def get_dimension(self):
        return Nimportant
    
    def change_basis(self,basis_func):
        pass