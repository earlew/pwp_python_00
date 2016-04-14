"""
This is script contains the helper functions for PWP.py
"""

import numpy as np
import seawater as sw
import matplotlib.pyplot as plt
from IPython.core.debugger import Tracer
debug_here = Tracer()

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
    
    #unpack arguments (TODO:not very elegant. Will fix in future)
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
    temp, sal, dens, uvel, vvel = remove_si(temp, sal, dens, uvel, vvel)
    
    #find ml index
    ml_thresh = 1e-4
    ml_idx = np.flatnonzero(np.diff(dens)>ml_thresh)[0] #finds the first index that exceed ML threshold
    ml_idx = ml_idx+1
    
    #check to ensure that ML is defined
    assert ml_idx.size is not 0, "Error: Mixed layer depth is undefined."
    
    #get surf MLD
    mld = z[ml_idx]
    
    #rotate u,v do wind input, rotate again, apply mixing
    ang = -f*dt/2
    uvel, vvel = rot(uvel, vvel, ang)
    du = (taux/(mld*dens[0]))*dt
    dv = (tauy/(mld*dens[0]))*dt
    uvel[:ml_idx] = uvel[:ml_idx]+du
    vvel[:ml_idx] = vvel[:ml_idx]+dv
    
    #apply drag to current 
    #Original comment: this is a horrible parameterization of inertial-internal wave dispersion
    if ucon > 1e-10:
        uvel = uvel*(1-dt*ucon) 
        vvel = vvel*(1-dt*ucon) 
    
    uvel, vvel = rot(uvel, vvel, ang)
    
    #Apply Bulk Richardson number instability form of mixing (as in PWP)
    if rb > 1e-5:
        temp, sal, dens, uvel, vvel = bulk_mix(t, s, d, u, v, dz, g, rg, zlen, ml_idx)
    
    if rc > 0:
        temp, sal, dens, uvel, vvel = grad_mix(temp, sal, dens, uvel, vvel, dz, g, rg, zlen)
    
    
    
def remove_si(temp, sal, dens, uvel, vvel):
    
    # Find and relieve static instability that may occur in the
    # density array d. This simulates free convection.
    # ml_index is the index of the depth of the surface mixed layer after adjustment,
    
    while stat_unstable:
        
        dens_diff = np.diff(dens)
        if np.any(dens_diff<0):
            stat_unstable=True
            first_inst_idx = np.flatnonzero(dens_diff<0)[0]
            initial_dens = dens
            (temp, sal, dens, uvel, vvel) = mix5(temp, sal, dens, uvel, vvel, first_inst_idx+1)
            
            #plot density
            plt.figure(num=86)
            plt.clf()
            plt.plot(initial_dens, 'b-')
            plt.plot(dens, 'r-')
            plt.show()
            
        else:
            stat_unstable = False
            
    return temp, sal, dens, uvel, vvel
    
    
def mix5(t, s, d, u, v, j):
    
    #This subroutine mixes the arrays t, s, u, v down to level j.
    j = j+1 #so that the j-th layer is included in the mixing
    t[:j] = np.mean(t[:j])
    s[:j] = np.mean(s[:j])
    d[:j] = sw.dens0(s[:j], t[:j])
    u[:j] = np.mean(u[:j])
    v[:j] = np.mean(v[:j])
    
    return t, s, d, u, v
    
        
def rot(u, v, ang):
    
    #This subroutine rotates the vector (u,v) through an angle, ang
    r = (u+1j*v)*np.exp(1j*ang)
    u = r.real
    v = r.imag
    
    return u, v   
    

def bulk_mix(t, s, d, u, v, dz, g, rg, zlen, ml_idx):
    #sub-routine to do bulk richardson mixing
    
    rvc = rb #critical rich number??
    
    for j in xrange(ml_ix, zlen+1):
    	h 	= z[j]
    	dd 	= (d[j]-d[1])/d[1];
    	dv 	= (u[j]-u[1])**2+(v[j]-v[1])**2;
    	if dv == 0:
    		rv = np.inf
    	else:
    		rv = g*h*dd/dv

    	if rv > rvc:
    		break
    	else:
    		t, s, d, u, v = mix5(t, s, d, u, v, j)
            
    return t, s, d, u, v

def grad_mix(t, s, d, u, v, dz, g, rg, zlen):
    
    #copied from source script:
    # %  This function performs the gradeint Richardson Number relaxation
    # %  by mixing adjacent cells just enough to bring them to a new
    # %  Richardson Number.      
    # %  Compute the gradeint Richardson Number, taking care to avoid dividing by
    # %  zero in the mixed layer.  The numerical values of the minimum allowable
    # %  density and velocity differences are entirely arbitrary, and should not
    # %  effect the calculations (except that on some occasions they evidently have!)
    
    rc = rg #critical rich. number
    j1 = 0
    j2 = zlen-1
    j_range = np.arange(ji,j2)
    
    while 1:
        #TODO: find a better way to do implement this loop
        
        r = np.zeros(len(j_range),)
        
        for j in j_range:
            # if j < 0:
            #     debug_here()
            
    		dd = (d[j+1]-d[j])/d[j]
    		dv = (u[j+1]-u[j])**2+(v[j+1]-v[j])**2
    		if dv == 0:
    			r[j] = np.inf
    		else:
                #compute grad. rich. number
    			r[j] = g*dz*dd/dv
                
                
        #find the smallest value of r in the profile
        r_min = np.min(r)
        j_min_idx = np.argmin(r)
        
        #Check to see whether the smallest r is critical or not.
        if r_min > rc:
            return t, s, d, u, v
            
        #Mix the cells js and js+1 that had the smallest Richardson Number
        t, s, d, u, v = stir(t, s, d, u, v, rc, j_min_idx)
        
        #recompute the rich number over the part of the profile that has changed
    	j1 = j_min_idx-2
    	if j1 < 1:
    		 j1 = 1
    	
    	j2 = j_min_idx+2
    	if j2 > zlen-1:
    		 j2 = zlen-1
             
    return t, s, d, u, v
        
        
def stir(t, s, d, u, v, rc, j):
    
    #copied from source script:
    
    # %  This subroutine mixes cells j and j+1 just enough so that
    # %  the Richardson number after the mixing is brought up to
    # %  the value rnew. In order to have this mixing process
    # %  converge, rnew must exceed the critical value of the
    # %  richardson number where mixing is presumed to start. If
    # %  r critical = rc = 0.25 (the nominal value), and r = 0.20, then
    # %  rnew = 0.3 would be reasonable. If r were smaller, then a
    # %  larger value of rnew - rc is used to hasten convergence.
    #
    # %  This subroutine was modified by JFP in Sep 93 to allow for an
    # %  aribtrary rc and to achieve faster convergence.
    
    #TODO: This needs better commenting
    rcon = 0.02+(rc-r)/2
    rnew = rc+rcon/5.
    f = 1-r/rnew
    
    dt = (t[j+1]-t[j])*f/2.
    t[j+1] = t[j+1]-dt
    t[j] = t[j]+dt
    
    ds = (s[j+1]-s[j])*f/2.
    s[j+1] = s[j+1]-ds
    s[j] = s[j]+ds
    d[j:j+1] = sw.dens0(s[j:j+1], t[j:j+1])
    
    du = (u[j+1]-u[j])*f/2
    u[j+1] = u[j+1]-du
    u[j] = u[j]+du
    
    dv = (v[j+1]-v[j])*f/2
    v[j+1] = v[j+1]-dv
    v[j] = v[j]+dv
    
    return t, s, d, u, v
    
    

            
            
    


















