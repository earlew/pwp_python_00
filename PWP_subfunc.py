"""
This is script contains the helper functions for PWP.py
"""

import numpy as np
import seawater as sw
import matplotlib.pyplot as plt

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

def pwpgo(forc, params, coords, temp, sal, uvel, vvel, dens):
    """
    Old doc string: this is where the model dynamics "happen"
    """
    
    #unpack arguments (not very elegant. Will fix in future)
    qi = forc['qi']
    q0 = forc['qo']
    emp = forc['emp']
    taux = forc['taux']
    tauy = forc['tauy']
    absrb = forc['absrb']
    
    z = coords['z']
    dz = coords['dz']
    dt = coords['dt']
    zlen = coords['zlen']
    
    rb = params['rb']
    rg = params['rg']
    f = params['f']
    cpw = params['cpw']
    g = params['g']
    ucon = params['ucon']
    
    ###Absorb solar radiation and FWF in surf layer ###
    
    #save initial T,S
    temp_old = temp[0]
    sal_old = sal[0]
    
    #update layer 1 temp and sal
    temp[0] = temp[0] + (qi*absrb[0]-qo)*dt/(dz*dens[0]*cpw)
    sal[0] = sal[0]/(1-emp*dt/dz)
    
    #check if temp is less than freezing point
    T_fz = sw.fp(sal_old, 1) #why 1? Need to recheck
    if temp[0] < T_fz:
        temp[0] = T_fz
        
    ###Absorb rad. at depth
    temp[1:] = temp[1:] + qi*absrb[1:]*dt/(dz*dens[1:]*cpw)
    
    #compute new density
    dens_new = sw.dens0(sal, temp)
    
    #relieve static instability
    
def remove_si(temp, sal, dens, uvel, vvel):
    
    # Find and relieve static instability that may occur in the
    # density array d. This simulates free convection.
    # ml_index is the index of the depth of the surface mixed layer after adjustment,
    
    while stat_unstable:
        
        dens_diff = np.diff(dens)
        if np.any(dens_diff<0):
            stat_unstable=True
            first_inst_ind = np.flatnonzero(dens_diff<0)[0]
            initial_dens = dens
            (temp, sal, dens, uvel, vvel) = mix5(temp, sal, dens, uvel, vvel, first_inst_ind+1)
            
            #plot density
            plt.figure(num=86)
            plt.clf()
            plt.plot(initial_dens, 'b-')
            plt.plot(dens, 'r-')
            plt.show()
            
        else:
            stat_unstable = False
            
        
def mix5(temp, sal, dens, uvel, vvel, j):
    
    #This subroutine mixes the arrays t, s, u, v down to level j.
    j = j+1 #so that the j-th layer is included in the mixing
    temp[:j] = np.mean(temp[:j])
    sal[:j] = np.mean(sal[:j])
    dens[:j] = sw.dens0(sal[:j], temp[:j])
    uvel[:j] = np.mean(uvel[:j])
    vvel[:j] = np.mean(vvel[:j])
    
    return temp, sal, dens, uvel, vvel
    
        
    
    




















