# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""

import numpy as np
import scipy.linalg as lin
import scipy.sparse as spar
from pyquantize.systematics.system import System
import matplotlib.pyplot as plt

class Coupler(System):
        
    def __init__(self,subsystems,name,coupling,**params):
        System.__init__(self,**params)
        
        self._name=name
        self._coupling=coupling
        
        self._ordered_ss=subsystems
        self._named_ss={}
        for name,ss in subsystems:
            self._named_ss.update({name:ss})
        self._params_updated()        
        
    def _params_updated(self):
        self._H_unc=0
        for name, sys in self._ordered_ss:
            if name in self._params:
                sys.set_params(**(self._params[name]))
            self._H_unc+= self.eval_op(name + '.H')
        H=self._H_unc+self.eval_op(self._coupling)
        self._energies,self._eigvecs=lin.eigh(H)     
            
    def _eval_simple_op(self,operator):
        if operator=='H':
            return np.matrix(np.diag(self._energies))
        else:
            sub_name,actual_op=operator.split('.')
    
            result=np.matrix([1])        
            for name,ss in self._ordered_ss:
                if name==sub_name:
                    result=np.kron(result,ss.eval_op(actual_op))
                else:
                    result=np.kron(result,ss.eval_op('eye'))
            return result
        return System._eval_simple_op(operator)
        
    def get_energies(self,n=None):
        return self._energies[:n]
        
    def get_dimension(self):
        return np.size(self._energies)
    
    def get_eigvecs(self):
        return self._eigvecs()
    
    def get_projectors(self,grouper):    
        
        # Get projectors onto pure grouper states
        N_g=self._named_ss[grouper].get_dimension()
        grouper_projectors=[]
        for n_g in range(N_g):
            grouper_projectors.append(
                self.eval_op(grouper + '.project' + str(n_g) ))
        return grouper_projectors
    
    def label_energies_by_grouper(self,grouper,which_energies='perturbed'):
        if which_energies=='both':
            return self.label_energies_by_grouper(grouper,'unperturbed'),\
                self.label_energies_by_grouper(grouper,'perturbed')
        projs=self.get_projectors(grouper)
        N=np.size(self._energies)
        
        def label(energies,states):
            levels=[]
            for n in range(N):
                state= states[:,n]
                energy=energies[n]
                mean_grouper= np.sum(map(\
                    lambda (n_g,proj): n_g*lin.norm(proj*state)**2,
                    enumerate(projs)))
                levels.append([mean_grouper,energy])
            levels.sort(key=lambda elt: elt[1])
            return levels        
        
        if which_energies=='unperturbed':
            return label(np.diag(self._H_unc),np.matrix(np.eye(N)))
        else:
            return label(self._energies,np.matrix(self._eigvecs))
    
    def plot_perturbation(self,grouper,Nlevels=None):
        
        unperturbed, perturbed=self.label_energies_by_grouper(grouper,'both')
        width=.25
        plt.figure()
        for g,e in unperturbed[:Nlevels]:
            plt.plot([g-width,g+width],[e,e],'b')
        for g,e in perturbed[:Nlevels]:
            plt.plot([g-width,g+width],[e,e],'--r')
        plt.show()

    def get_grouped_spectrum(self,grouper,which_energies='perturbed'):
        if which_energies=='both':
            return self.get_grouped_spectrum(grouper,'unperturbed'),\
                self.get_grouped_spectrum(grouper,'perturbed')
        
        specs=[]
        for n in range(self._named_ss[grouper].get_dimension()):
            spec=[e for g,e in\
                self.label_energies_by_grouper(grouper,which_energies)      
                if abs(n-g)<.5]
            specs.append(spec)
        return specs
        
plt.close('all')
from pyquantize.systematics.system_sho import System_SHO
f1=5;  mysho1=System_SHO(freq=f1,Nbasis=2)
f2=5.1;mysho2=System_SHO(freq=f2,Nbasis=40)
myc=Coupler([['one',mysho1],['two',mysho2]],'couple',\
    [lambda a1,ad1,a2,ad2: .01*(a1+ad1)*(a2+ad2),
         ['one.a','one.a^T','two.a','two.a^T']])
myc.eval_op('one.a')
myc.eval_op('one.project1')
myc.plot_perturbation('one',10)
grouped=myc.get_grouped_spectrum('one')
print (grouped[0][1]-grouped[0][0])-f2
print (grouped[1][1]-grouped[1][0])-f2