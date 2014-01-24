# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:15:09 2014

@author: sam
"""
import numpy as np
import re

class System(object):
        
    def __init__(self,**params):
        self._params=params
        
    def set_params(self,**params):
        self._params.update(params)
        self._params_updated()
        
    def _params_updated(self):
        pass
    
    def eval_op(self,operator):
            
        if hasattr(operator[0],'__call__'):
            ops=map(lambda op: self.eval_op(op),operator[1])
            return operator[0](*ops)
        elif operator=='eye':
            return np.matrix(np.eye(self.get_dimension()))
            
        m=re.compile('project(\d+)').match(operator)
        if m is not None:
            state=int(m.group(1))
            return np.diag((np.arange(self.get_dimension())==state)*1)
        
        return self._eval_simple_op(operator)
            
    def _eval_simple_op(self,operator):
        print "OPERATOR NOT RECOGNIZED"
        
    def get_energies(self):
        pass
        
    def get_eigvecs(self):
        pass
        
    def get_dimension(self):
        pass