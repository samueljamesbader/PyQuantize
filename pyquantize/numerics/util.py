# -*- coding: utf-8 -*-
"""
Module: qsim.numerics.util
Author: Sam Bader

Various utility functions for numerics
"""
import numpy as np
import matplotlib.pyplot as plt


def force_array(mat):
    ''' Forces a matrix into a 1-D array.'''
    return np.squeeze(np.array(mat))
    
def force_plot(x,y,*args,**kwargs):
    plt.plot(force_array(x),force_array(y),*args,**kwargs)
    plt.show()