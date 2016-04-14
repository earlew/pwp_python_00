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
    taux = met_dset_intp['tx'] #x-component of wind stress
    tauy = met_dset_intp['ty'] #y-component of wind stress
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
        
    #get profile variables
    temp = prof_dset_intp['t'] #profile temperature
    sal = prof_dset_intp['s'] #profile salinity
    dens = sw.dens0(sal, temp)
     
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
    uvel = np.zeros((zlen, tlen)) #east velocity m/s
    vvel = np.zeros((zlen, tlen)) #north velocity m/s
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
    
    #pre-allocate loop vbles
    sal_n = np.zeros((zlen, tlen))
    temp_n = np.zeros((zlen, tlen))
    vvel_n = np.zeros((zlen, tlen))
    uvel_n = np.zeros((zlen, tlen))
    mld_n = np.zeros((tlen,))
    
    for n in xrange(1,tlen):
        print 'Loop iter. %s' %n
        
        #package input variables (TODO: find a better way to do this. Maybe use classes?)
        nm1 = n-1
        forc = {}
        forc['q_in'] = q_in[nm1]
        forc['q_out'] = q_out[nm1]
        forc['emp'] = emp[nm1]
        forc['taux'] = taux[nm1]
        forc['tauy'] = tauy[nm1]
        forc['absrb'] = absrb
        
        coords['z'] = z
        coords['dz'] = dz
        coords['dt'] = dt
        coords['zlen'] = zlen
        
        params['rb'] = rb
        params['rg'] = rg
        params['f'] = f
        params['cpw'] = cpw
        params['g'] = g
        params['ucon'] = ucon
        
        sal_n[:,n], temp_n[:,n], uvel_n[:,n], vvel_n[:,n], mld_n[n] = sbf.pwpgo(forc, params, coords, temp, sal, uvel, vvel, dens, n)

        #apply vertical diffusion
        if rkz > 0:
            #diffusion
            #this code block is incomplete in the source script
            
        #Diagnostic plots
        if diagnostics == 1:
            
            #plot depth int. KE and momentum
            plt.figure(num=1)
            
            plt.subplot(211)
            plt.plot(time[n]-time[0], np.trapz(z, 0.5*d*(uvel[:,n]**2+vvel[:,n]**2)), 'b.')
            plt.grid(True)
            if n==1:
                plt.title('Depth integrated KE')
            
            plt.subplot(212)
            plt.plot(time[n]-time[0], np.trapz(z, d*np.sqrt(uvel[:,n]**2+vvel[:,n]**2)), 'b.')
            plt.grid(True)
            if n==1:
                plt.title('Depth integrated Mom.')
                #plt.get_current_fig_manager().window.wm_geometry("400x600+20+40")
                
            #plot T,S and U,V
            ax1 = plt.subplot2grid((4,1), (0, 0), colspan=2)
            ax1.plot(uvel[:,n], z, 'b', label='uvel')
            ax1.plot(vvel[:,n], z, 'r', label='vvel')
            ax1.grid(True)
            ax1.legend()
            
            ax2 = plt.subplot2grid((4,1), (0, 2), colspan=1)
            ax2.plot(temp[:,n], z, 'b')
            ax2.grid(True)
            ax2.set_xlabel('Temp.')
            
            ax3 = plt.subplot2grid((4,1), (0, 3), colspan=1)
            ax2.plot(sal[:,n], z, 'b')
            ax3.grid(True)
            ax3.set_xlabel('Salinity.')
            
            
            
            
        
if __name__ == "__main__":
    
    main()    
    
    
    
    
    
    
    





