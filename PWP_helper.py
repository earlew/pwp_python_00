"""
This module contains the helper functions to assist with the running and analysis of the
PWP model.
"""

import numpy as np
import seawater as sw
import matplotlib.pyplot as plt
import PWP
from datetime import datetime
import warnings

#warnings.filterwarnings("error")
#warnings.simplefilter('error', RuntimeWarning)

#from IPython.core.debugger import set_trace
#debug_here = set_trace

def run_demo1():
    """
    Example script of how to run the PWP model. 
    This run uses summertime data from the Beaufort gyre
    """
    
    #ds = xr.Dataset({'t': (['z'], prof_dset['t'].values), 's': (['z'], prof_dset['s'].values), 'd': (['z'], prof_dset['d'].values), 'z': (['z'], prof_dset['z'].values), 'lat': 74.0})
    
    forcing_fname = 'beaufort_met.nc'
    prof_fname = 'beaufort_profile.nc'
    print("Running Test Case 1 with data from Beaufort gyre...")
    forcing, pwp_out = PWP.run(met_data=forcing_fname, prof_data=prof_fname, suffix='demo1_nodiff', save_plots=True, diagnostics=False)
    
    #debug_here()
def run_demo2(winds_ON=True, emp_ON=True, heat_ON=True, drag_ON=True):
    
    """
    Example script of how to run the PWP model.
    This run uses summertime data from the Atlantic sector of the Southern Ocean
    """
    
    forcing_fname = 'SO_met_30day.nc'
    prof_fname = 'SO_profile1.nc'
    print("Running Test Case 2 with data from Southern Ocean...")
    p={}
    p['rkz']=1e-6
    p['dz'] = 2.0 
    p['max_depth'] = 500.0 
    p['rg'] = 0.25 # to turn off gradient richardson number mixing set to 0. (code runs much faster)
    p['winds_ON'] = winds_ON
    p['emp_ON'] = emp_ON
    p['heat_ON'] = heat_ON
    p['drag_ON'] = drag_ON
    
    if emp_ON: 
        emp_flag=''
    else:
        emp_flag='_empOFF'
        
    if winds_ON: 
        winds_flag=''
    else:
        winds_flag='_windsOFF'
        
    if heat_ON: 
        heat_flag=''
    else:
        heat_flag='_heatingOFF'
        
    if drag_ON:
        drag_flag=''
    else:
        drag_flag='_dragOFF'
        
    suffix = 'demo2_1e6diff%s%s%s%s' %(winds_flag, emp_flag, heat_flag, drag_flag)
    forcing, pwp_out = PWP.run(met_data=forcing_fname, prof_data=prof_fname, suffix=suffix, save_plots=True, param_kwds=p)
     

def set_params(lat, dt=3., dz=1., max_depth=100., mld_thresh=1e-4, dt_save=1., rb=0.65, rg=0.25, rkz=0., beta1=0.6, beta2=20.0, heat_ON=True, winds_ON=True, emp_ON=True, drag_ON=True):
    
    """
    This function sets the main paramaters/constants used in the model.
    These values are packaged into a dictionary, which is returned as output.
    Definitions are listed below.
    
    CONTROLS (default values are in [ ]):
    lat: latitude of profile
    dt: time-step increment. Input value in units of hours, but this is immediately converted to seconds.[3 hours]
    dz: depth increment (meters). [1m]
    max_depth: Max depth of vertical coordinate (meters). [100]
    mld_thresh: Density criterion for MLD (kg/m3). [1e-4] 
    dt_save: time-step increment for saving to file (multiples of dt). [1]
    rb: critical bulk richardson number. [0.65]
    rg: critical gradient richardson number. [0.25]
    rkz: background vertical diffusion (m**2/s). [0.]
    beta1: longwave extinction coefficient (meters). [0.6] 
    beta2: shortwave extinction coefficient (meters). [20] 
    winds_ON: True/False flag to turn ON/OFF wind forcing. [True]
    emp_ON: True/False flag to turn ON/OFF freshwater forcing. [True]
    heat_ON: True/False flag to turn ON/OFF surface heat flux forcing. [True]
    drag_ON: True/False flag to turn ON/OFF current drag due to internal-inertial wave breaking. [True]
    
    OUTPUT is dict with fields containing the above variables plus the following:
    dt_d: time increment (dt) in units of days
    g: acceleration due to gravity [9.8 m/s^2]
    cpw: specific heat of water [4183.3 J/kgC]
    f: coriolis term (rad/s). [sw.f(lat)]
    ucon: coefficient of inertial-internal wave dissipation (s^-1) [0.1*np.abs(f)]
    """
    params = {}
    params['dt'] = 3600.0*dt
    params['dt_d'] = params['dt']/86400.
    params['dz'] = dz
    params['dt_save'] = dt_save
    params['lat'] = lat
    params['rb'] = rb
    params['rg'] = rg
    params['rkz'] = rkz
    params['beta1'] = beta1
    params['beta2'] = beta2
    params['max_depth'] = max_depth
  
    params['g'] = 9.81
    params['f'] = sw.f(lat)
    params['cpw'] = 4183.3
    params['ucon'] = (0.1*np.abs(params['f']))
    params['mld_thresh'] = mld_thresh
    
    params['winds_ON'] = winds_ON
    params['emp_ON'] = emp_ON
    params['heat_ON'] = heat_ON
    params['drag_ON'] = drag_ON
    
    return params
    
    

def prep_data(met_dset, prof_dset, params):
    
    """
    This function prepares the forcing and profile data for the model run.
    
    Below, the surface forcing and profile data are interpolated to the user defined time steps
    and vertical resolutions, respectively. Secondary quantities are also computed and packaged 
    into dictionaries. The code also checks that the time and vertical increments meet the 
    necessary stability requirements.
    
    Lastly, this function initializes the numpy arrays to collect the model's output.
    
    INPUT:
    met_data: dictionary-like object with forcing data. Fields should include: 
            ['time', 'sw', 'lw', 'qlat', 'qsens', 'tx', 'ty', 'precip']. These fields should 
            store 1-D time series of the same length. 
            
            The model expects positive heat flux values to represent ocean warming. The time
            data field should contain a 1-D array representing fraction of day. For example, 
            for 6 hourly data, met_data['time'] should contain a number series that increases
            in steps of 0.25, such as np.array([1.0, 1.25, 1.75, 2.0, 2.25...]).

            See https://github.com/earlew/pwp_python#input-data for more info about the
            expect intput data. 
    
            TODO: Modify code to accept met_data['time'] as an array of datetime objects
    
            
    prof_data: dictionary-like object with initial profile data. Fields should include:
            ['z', 't', 's', 'lat']. These represent 1-D vertical profiles of temperature,
            salinity and density. 'lat' is expected to be a length=1 array-like object. e.g. 
            prof_data['lat'] = [25.0]
            
    params: dictionary-like object with fields defined by set_params function
    
    OUTPUT:
    
    forcing: dictionary with interpolated surface forcing data. 
    pwp_out: dictionary with initialized variables to collect model output.
    """
    
    #create new time vector with time step dt_d
    #time_vec = np.arange(met_dset['time'][0], met_dset['time'][-1]+params['dt_d'], params['dt_d']) 
    time_vec = np.arange(met_dset['time'][0], met_dset['time'][-1], params['dt_d']) 
    tlen = len(time_vec)
    
    #debug_here()
    
    #interpolate surface forcing data to new time vector
    from scipy.interpolate import interp1d
    forcing = {} 
    for vname in met_dset:
        p_intp = interp1d(met_dset['time'], met_dset[vname], axis=0)
        forcing[vname] = p_intp(time_vec)
        
    #interpolate E-P to dt resolution (not sure why this has to be done separately)
    evap_intp = interp1d(met_dset['time'], met_dset['qlat'], axis=0, kind='nearest', bounds_error=False)
    evap = (0.03456/(86400*1000))*evap_intp(np.floor(time_vec)) #(meters per second?)
    emp = np.abs(evap) - np.abs(forcing['precip'])
    emp[np.isnan(emp)] = 0.
    forcing['emp'] = emp 
    forcing['evap'] = evap 
    
    if params['emp_ON'] == False:
        print("WARNING: E-P is turned OFF.")
        forcing['emp'][:] = 0.0
        forcing['precip'][:] = 0.0
        forcing['evap'][:] = 0.0
        
    if params['heat_ON'] == False:
        print("WARNING: Surface heating is turned OFF.")
        forcing['sw'][:] = 0.0
        forcing['lw'][:] = 0.0
        forcing['qlat'][:] = 0.0
        forcing['qsens'][:] = 0.0
          
    
    #define q_in and q_out (positive values should mean ocean warming)
    forcing['q_in'] = forcing['sw'] #heat flux into ocean
    forcing['q_out'] = -(forcing['lw'] + forcing['qlat'] + forcing['qsens']) 
    
    #add time_vec to forcing
    forcing['time'] = time_vec
    
    if params['winds_ON'] == False:
        print("Winds are set to OFF.")
        forcing['tx'][:] = 0.0
        forcing['ty'][:] = 0.0
           
    #define depth coordinate, but first check to see if profile max depth
    #is greater than user defined max depth
    zmax = max(prof_dset.z)
    if zmax < params['max_depth']:
        depth = zmax
        print('Profile input shorter than depth selected, truncating to %sm' %depth)
        
    
    #define new z-coordinates
    init_prof = {}
    init_prof['z'] = np.arange(0, params['max_depth']+params['dz'], params['dz'])
    zlen = len(init_prof['z'])
    
    #compute absorption and incoming radiation (function defined in PWP_model.py)
    absrb = PWP.absorb(params['beta1'], params['beta2'], zlen, params['dz']) #(units unclear)
    dstab = params['dt']*params['rkz']/params['dz']**2 #courant number  
    if dstab > 0.5:
        print("WARNING: unstable CFL condition for diffusion! dt*rkz/dz**2 > 0.5.")
        print("To fix this, try to reduce the time step or increase the depth increment.")
        inpt = eval(input("Proceed with simulation? Enter 'y'or 'n'. "))
        if inpt == 'n':
            raise ValueError("Please restart PWP.m with a larger dz and/or smaller dt. Exiting...")
        
    forcing['absrb'] = absrb
    params['dstab'] = dstab
    
    #check depth resolution of profile data
    prof_incr = np.diff(prof_dset['z']).mean()
    # if params['dz'] < prof_incr/5.:
    #     message = "Specified depth increment (%s m), is much smaller than mean profile resolution (%s m)." %(params['dz'], prof_incr)
    #     warnings.warn(message)
        
    
    #debug_here()
    #interpolate profile data to new z-coordinate
    from scipy.interpolate import InterpolatedUnivariateSpline 
    for vname in prof_dset:
        if vname == 'lat' or vname=='lon':
            continue
        else:
            #first strip nans
            not_nan = np.logical_not(np.isnan(prof_dset[vname]))
            indices = np.arange(len(prof_dset[vname]))
            #p_intp = interp1d(prof_dset['z'], prof_dset[vname], axis=0, kind='linear', bounds_error=False)
            #interp1d doesn't work here because it doesn't extrapolate. Can't have Nans in interpolated profile
            p_intp = InterpolatedUnivariateSpline(prof_dset['z'][not_nan], prof_dset[vname][not_nan], k=1)
            init_prof[vname] = p_intp(init_prof['z'])    
        
    #get profile variables
    temp0 = init_prof['t'] #initial profile temperature
    sal0 = init_prof['s'] #intial profile salinity
    dens0 = sw.dens0(sal0, temp0) #intial profile density       
    
    #initialize variables for output
    pwp_out = {}
    pwp_out['time'] = time_vec
    pwp_out['dt'] = params['dt']
    pwp_out['dz'] = params['dz']
    pwp_out['lat'] = params['lat']
    pwp_out['z'] = init_prof['z']
    
    tlen = int(np.floor(tlen/params['dt_save']))
    arr_sz = (zlen, tlen)
    pwp_out['temp'] = np.zeros(arr_sz)
    pwp_out['sal'] = np.zeros(arr_sz)
    pwp_out['dens'] = np.zeros(arr_sz)
    pwp_out['uvel'] = np.zeros(arr_sz)
    pwp_out['vvel'] = np.zeros(arr_sz)
    pwp_out['mld'] = np.zeros((tlen,))
    
    #use temp, sal and dens profile data for the first time step
    pwp_out['sal'][:,0] = sal0
    pwp_out['temp'][:,0] = temp0
    pwp_out['dens'][:,0] = dens0
    
    return forcing, pwp_out, params
    
def livePlots(pwp_out, n):
    
    """
    function to make live plots of the model output.
    """
    
    #too lazy to re-write the plotting code, so i'm just going to unpack pwp_out here:
    time = pwp_out['time']
    uvel = pwp_out['uvel']
    vvel = pwp_out['vvel']
    temp = pwp_out['temp']
    sal = pwp_out['sal']
    dens = pwp_out['dens']
    z = pwp_out['z']


    #plot depth int. KE and momentum
    plt.figure(num=1)

    plt.subplot(211)
    plt.plot(time[n]-time[0], np.trapz(0.5*dens[:,n]*(uvel[:,n]**2+vvel[:,n]**2)), 'b.')
    plt.grid(True)
    if n==1:
        plt.title('Depth integrated KE')

    plt.subplot(212)
    plt.plot(time[n]-time[0], np.trapz(dens[:,n]*np.sqrt(uvel[:,n]**2+vvel[:,n]**2)), 'b.')
    plt.grid(True)
    plt.pause(0.05)
    plt.subplots_adjust(hspace=0.35)

    #debug_here()
    if n==1:
        plt.title('Depth integrated Mom.')
        #plt.get_current_fig_manager().window.wm_geometry("400x600+20+40")
    
    #plot T,S and U,V
    plt.figure(num=2, figsize=(12,6))
    ax1 = plt.subplot2grid((1,4), (0, 0), colspan=2)
    ax1.plot(uvel[:,n], z, 'b', label='uvel')
    ax1.plot(vvel[:,n], z, 'r', label='vvel')
    ax1.invert_yaxis()
    ax1.grid(True)
    ax1.legend(loc=3)    

    ax2 = plt.subplot2grid((1,4), (0, 2), colspan=1)
    ax2.plot(temp[:,n], z, 'b')
    ax2.grid(True)
    ax2.set_xlabel('Temp.')
    ax2.invert_yaxis()
    xlims = ax2.get_xlim()
    xticks = np.round(np.linspace(xlims[0], xlims[1], 4), 1)
    ax2.set_xticks(xticks)

    ax3 = plt.subplot2grid((1,4), (0, 3), colspan=1)
    ax3.plot(sal[:,n], z, 'b')
    ax3.set_xlabel('Salinity')
    ax3.grid(True)
    ax3.invert_yaxis()
    xlims = ax3.get_xlim()
    xticks = np.round(np.linspace(xlims[0], xlims[1], 4), 1)
    ax3.set_xticks(xticks)

    plt.pause(0.05)

    plt.show()

def makeSomePlots(forcing, pwp_out, time_vec=None, save_plots=False, suffix=''):
    
    """
    TODO: add doc file
    Function to make plots of the results once the model iterations are complete.
    
    """
    
    if len(suffix)>0 and suffix[0] != '_':
            suffix = '_%s' %suffix
    
    #plot summary of ML evolution
    fig, axes = plt.subplots(3,1, sharex=True, figsize=(7.5,9))
    
    if time_vec is None:
        tvec = pwp_out['time']
    else:
        tvec = time_vec
    
    axes = axes.flatten()
    ##plot surface heat flux
    axes[0].plot(tvec, forcing['lw'], label='$Q_{lw}$')
    axes[0].plot(tvec, forcing['qlat'], label='$Q_{lat}$')
    axes[0].plot(tvec, forcing['qsens'], label='$Q_{sens}$')
    axes[0].plot(tvec, forcing['sw'], label='$Q_{sw}$')
    axes[0].hlines(0, tvec[0], pwp_out['time'][-1], linestyle='-', color='0.3')
    axes[0].plot(tvec, forcing['q_in']-forcing['q_out'], ls='-', lw=2, color='k', label='$Q_{net}$')   
    axes[0].set_ylabel('Heat flux (W/m2)')
    axes[0].set_title('Heat flux into ocean')
    axes[0].grid(True)
    #axes[0].set_ylim(-500,300)
    
    axes[0].legend(loc=0, ncol=2, fontsize='smaller')
    
    
    ##plot wind stress
    axes[1].plot(tvec, forcing['tx'], label=r'$\tau_x$')
    axes[1].plot(tvec, forcing['ty'], label=r'$\tau_y$')
    axes[1].hlines(0, tvec[0], pwp_out['time'][-1], linestyle='--', color='0.3')
    axes[1].set_ylabel('Wind stress (N/m2)')
    axes[1].set_title('Wind stress')
    axes[1].grid(True)
    axes[1].legend(loc=0, fontsize='medium')
    
    
    ## plot freshwater forcing
    # emp_mmpd = forcing['emp']*1000*3600*24 #convert to mm per day
    # axes[2].plot(tvec, emp_mmpd, label='E-P')
    # axes[2].hlines(0, tvec[0], pwp_out['time'][-1], linestyle='--', color='0.3')
    # axes[2].set_ylabel('Freshwater forcing (mm/day)')
    # axes[2].set_title('Freshwater forcing')
    # axes[2].grid(True)
    # axes[2].legend(loc=0, fontsize='medium')
    # axes[2].set_xlabel('Time (days)')
    
    emp_mmpd = forcing['emp']*1000*3600*24 #convert to mm per day
    evap_mmpd = forcing['evap']*1000*3600*24 #convert to mm per day
    precip_mmpd = forcing['precip']*1000*3600*24 #convert to mm per day
    axes[2].plot(tvec, precip_mmpd, label='$P$', lw=1, color='b')
    axes[2].plot(tvec, evap_mmpd, label='$-E$', lw=1, color='r')
    axes[2].plot(tvec, emp_mmpd, label='$|E| - P$', lw=2, color='k')
    axes[2].hlines(0, tvec[0], tvec[-1], linestyle='--', color='0.3')
    axes[2].set_ylabel('Freshwater forcing (mm/day)')
    axes[2].set_title('Freshwater forcing')
    axes[2].grid(True)
    axes[2].legend(loc=1, fontsize=8, ncol=2)
    axes[2].set_xlabel('Time (days)')
    
    if save_plots:     
        plt.savefig('plots/surface_forcing%s.png' %suffix, bbox_inches='tight')
    
    
    ##plot temp and sal change over time
    fig, axes = plt.subplots(2,1, sharex=True)
    vble = ['temp', 'sal']
    units = ['$^{\circ}$C', 'PSU']
    #cmap = custom_div_cmap(numcolors=17)
    cmap = plt.cm.rainbow
    for i in range(2):
        ax = axes[i]
        im = ax.contourf(pwp_out['time'], pwp_out['z'], pwp_out[vble[i]], 15, cmap=cmap, extend='both')
        ax.set_ylabel('Depth (m)')
        ax.set_title('Evolution of ocean %s (%s)' %(vble[i], units[i]))
        ax.invert_yaxis()   
        cb = plt.colorbar(im, ax=ax, format='%.1f')
     
    ax.set_xlabel('Days')   
    
    
    ## plot initial and final T-S profiles
    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    
    plt.figure()
    host = host_subplot(111, axes_class=AA.Axes)
    host.invert_yaxis()
    par1 = host.twiny() #par for parasite axis
    host.set_ylabel("Depth (m)")
    host.set_xlabel("Temperature ($^{\circ}$C)")
    par1.set_xlabel("Salinity (PSU)")
    
    p1, = host.plot(pwp_out['temp'][:,0], pwp_out['z'], '--r', label='$T_i$')
    host.plot(pwp_out['temp'][:,-1], pwp_out['z'], '-r', label='$T_f$')
    p2, = par1.plot(pwp_out['sal'][:,0], pwp_out['z'], '--b', label='$S_i$')
    par1.plot(pwp_out['sal'][:,-1], pwp_out['z'], '-b', label='$S_f$')
    host.grid(True)
    host.legend(loc=0, ncol=2)
    #par1.legend(loc=3)
    
    host.axis["bottom"].label.set_color(p1.get_color())
    host.axis["bottom"].major_ticklabels.set_color(p1.get_color())
    host.axis["bottom"].major_ticks.set_color(p1.get_color())

    par1.axis["top"].label.set_color(p2.get_color())
    par1.axis["top"].major_ticklabels.set_color(p2.get_color())
    par1.axis["top"].major_ticks.set_color(p2.get_color())
    
    if save_plots:     
        plt.savefig('plots/initial_final_TS_profiles%s.png' %suffix, bbox_inches='tight')
        
    
    
    plt.show()
    

def custom_div_cmap(numcolors=11, name='custom_div_cmap', mincol='blue', midcol='white', maxcol='red'):
                    
    """ Create a custom diverging colormap with three colors
    
    Default is blue to white to red with 11 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    """

    from matplotlib.colors import LinearSegmentedColormap 
    
    cmap = LinearSegmentedColormap.from_list(name=name, 
                                             colors =[mincol, midcol, maxcol],
                                             N=numcolors)
    return cmap   