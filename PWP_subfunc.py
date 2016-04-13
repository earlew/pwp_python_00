"""
This is script contains the helper functions for PWP.py
"""

import numpy as np
import seawater as sw

def absorb(beta1, beta2, zlen, dz):
    
    # Compute solar radiation absorption profile. This
    # subroutine assumes two wavelengths, and a double
    # exponential depth dependence for absorption.
    # 
    # Subscript 1 is for red, non-penetrating light, and
    # 2 is for blue, penetrating light. rs1 is the fraction
    # assumed to be red.
    
    rs1 = 0.6
    rs2 = 1.0-rs1
    z1 = np.arange(0,zlen)*dz
    z2 = z1 + dz
    z1b1 = z1/beta1
    z2b1 = z2/beta1
    z1b2 = z1/beta2
    z2b2 = z2/beta2
    absrb = rs1*(np.exp(-z1b1)-np.exp(-z2b1))+rs2*(np.exp(-z1b2)-np.exp(-z2b2))
    
    return absrb

def pwpgo():
    
    pass



