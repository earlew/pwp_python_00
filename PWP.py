"""
Python implementation of the Price Weller Pinkel (PWP) mixed layer model.
This is based on the MATLAB implementation of the PWP model written by
Byron Kilbourne (University of Washington)
"""

#import universal modules
import numpy as np
import seawater as sw
import xray
import PWP_subfunc as sbf
from IPython.core.debugger import Tracer
reload(sbf)

debug_here = Tracer()

#start of main  script
def main():
    
    ##############################################
    #define model parameters
    ##############################################
    
    dt = 3600.0*3 #time-step increment (seconds)
    secsInDay = 86400.
    dt_d = dt/secsInDay
    dz = 1.  #depth increment (meters)
    max_depth = 100. #the depth to run (max depth grid)
    dt_save = 1. #time-step increment for saving to file (multiples of dt)
    lat = 74. #latitude (degrees)
    g = 9.8 #gravity (9.8 m/s^2)
    cpw	= 4183.3 #specific heat of water (4183.3 J/kgC)
    rb = 0.65 #critical bulk richardson number (0.65)
    rg = 0.25 #critical gradient richardson number (0.25)
    rkz	= 0. #background vertical diffusion (0) m^2/s
    beta1 = 0.6 #longwave extinction coefficient (0.6 m)
    beta2 = 20.  #shortwave extinction coefficient (20 m)
    f = sw.f(lat) #coriolis term (rad/s)
    ucon = (.1*np.abs(f)); #coefficient of inertial-internal wave dissipation (0) s^-1
    
    ##############################################
    #Prepare surface forcing and profile data
    ##############################################  
    
    met_dset = xray.open_dataset('met.nc')
    prof_dset = xray.open_dataset('profile.nc')

    #create new time vector with time step dt_d
    time_vec = np.arange(met_dset['time'][0], met_dset['time'][-1]+dt_d, dt_d) 
    tlen = len(time_vec)
    
    #interpolate surface forcing data to new time vector
    from scipy.interpolate import interp1d
    met_dset_intp = {} 
    for vname in met_dset:
        p_intp = interp1d(met_dset['time'], met_dset[vname], axis=0)
        met_dset_intp[vname] = p_intp(time_vec)
    
    #unpack interpolated forcing variables
    q_in = met_dset_intp['sw'] #heat flux into ocean
    q_out = met_dset_intp['lw'] + met_dset_intp['qlat'] + met_dset_intp['qsens'] #heat flux out of ocean
    tau_x = met_dset_intp['tx'] #x-component of wind stress
    tau_y = met_dset_intp['ty'] #y-component of wind stress
    precip = met_dset_intp['precip'] #precip 
    
    #define depth coordinate, but first check to see if profile max depth
    #is greater than user defined max depth
    zmax = max(prof_dset.z)
    if zmax < max_depth:
        depth = zmax
        print 'Profile input shorter than depth selected, truncating to %sm' %depth
    
    z = np.arange(0, max_depth+dz, dz)
    zlen = len(z)
    
    #check depth resolution of profile data
    prof_incr = np.diff(prof_dset['z']).mean()
    if dz < prof_incr/5.:
        inpt = input('Depth increment, dz, is much smaller than profile resolution. Is this okay? (y/n)')
        if inpt is 'n':
            print "Please restart PWP.m with a new dz >= %s" %prof_incr/5.
            return
            
    #interpolate profile data to new z-coordinate
    prof_dset_intp = {} 
    for vname in prof_dset:
        p_intp = interp1d(prof_dset['z'], prof_dset[vname], axis=0, kind='nearest', bounds_error=False)
        prof_dset_intp[vname] = p_intp(z)
        
    prof_temp = prof_dset_intp['t'] #profile temperature
    prof_sal = prof_dset_intp['s'] #profile salinity
    prof_dens = sw.dens0(prof_sal, prof_temp)
     
    #interpolate E-P to dt resolution
    evap = (0.03456/(86400*1000))*met_dset_intp['qlat'] #(meters?)
    emp = evap-precip
    
    #compute absorption and courant number
    absrb = sbf.absorb(beta1, beta2, zlen, dz) #absorption of icncoming radiation (units unclear)
    dstab = dt*rkz/dz**2 #courant number  
    if dstab > 0.5:
        print 'WARNING: !unstable CFL condition for diffusion!'
    
    ##############################################
    #Prepare variables for output 
    ##############################################
    u = np.zeros((zlen, tlen)) #east velocity m/s
    v = np.zeros((zlen, tlen)) #north velocity m/s
    mld = np.zeros((tlen,))
    
    #preallocate output dict
    pwp_output = {}
    pwp_output['dt'] = dt
    pwp_output['dz'] = dz
    pwp_output['lat'] = lat
    pwp_output['z'] = z
    pwp_output['time'] = np.zeros((zlen, np.floor(tlen/dt_save)))
    pwp_output['temp'] = np.zeros((zlen, np.floor(tlen/dt_save)))
    pwp_output['sal'] = np.zeros((zlen, np.floor(tlen/dt_save)))
    pwp_output['dens'] = np.zeros((zlen, np.floor(tlen/dt_save)))
    pwp_output['uvel'] = np.zeros((zlen, np.floor(tlen/dt_save)))
    pwp_output['vvel'] = np.zeros((zlen, np.floor(tlen/dt_save)))
    pwp_output['mld'] = np.zeros((zlen, np.floor(tlen/dt_save)))
    
    ##############################################
    #MODEL LOOP START 
    ##############################################
    
    for n in xrange(1,tlen):
        print 'Loop iter. %s' %n
        
if __name__ == "__main__":
    
    main()    
    
    
    
    
    
    
    





