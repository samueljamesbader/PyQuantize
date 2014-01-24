# -*- coding: utf-8 -*-
"""
Module: qsim.numerics.bases.util
Author: Convenient functions for examining bases

"""
import matplotlib.pyplot as plt

def check_basis(basis_mat):
    ''' Plots the norm of the basis vectors versus basis vector number '''

    norm=[(basis_mat[:,i].T*basis_mat[:,i])[0,0]\
        for i in range(basis_mat.shape[1])]
    plt.figure()
    plt.plot(norm)
    plt.axis([0,50,.9,1.1])