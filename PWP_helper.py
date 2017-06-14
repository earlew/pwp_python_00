"""
This module contains the helper functions to assist with the running and analysis of the
PWP model.
"""

import numpy as np
import seawater as sw
import matplotlib.pyplot as plt
import PWP
from datetime import datetime
import xarray as xr
import timeit
import pickle

import warnings
from IPython.core.debugger import Tracer
debug_here = Tracer()

def demo1():
    """
    Example script of how to run the PWP model. 
    This run uses summertime data from the Beaufort gyre
    """
    
    forcing_fname = 'beaufort_met.nc'
    prof_fname = 'beaufort_profile.nc' 
    print("Running Test Case 1 with data from Beaufort gyre...")
    forcing, pwp_out = run_PWP(met_data=forcing_fname, prof_data=prof_fname, makeLivePlots=False, suffix='demo1_nodiff', save_plots=True)
    
    return forcing, pwp_out

def demo2(dopt='pdens'):
    
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
    p['dopt'] = dopt
    warnings.simplefilter('error', UserWarning)
    forcing, pwp_out = run_PWP(met_data=forcing_fname, prof_data=prof_fname, suffix='demo2_1e6diff', save_plots=True, param_kwds=p)
    
    return forcing, pwp_out
    

def demo3():
    
    """
    Example script of how to run the PWP model.
    This run is initialized with an early winter profile from Maud Rise and uses NCEP fluxes. 
    """
    
    p={}
    p['rkz']=1e-6
    p['dz'] = 1
    p['dt'] = 1.5
    p['max_depth'] = 500
    p['ice_ON'] = True
    p['winds_ON'] = True
    p['emp_ON'] = False
    p['alpha'] = 0.95
    p['dopt'] = 'pdens'
    p['fix_alpha'] = True
    p['mld_thresh'] = 0.01
    p['use_Bulk_Formula'] = True
    p['qnet_offset'] = 0 #W/m2
    p['iceMod'] = 1 #1: use ice_model_0(), 1: ice_model_T()

    
    if p['ice_ON']:
        ice_str = '' 
    else:
        ice_str = '_noICE'
        
    if p['winds_ON']:
        wind_str = ''
    else:
        wind_str = '_noWINDS'
        
    if p['emp_ON']:
        emp_str = ''
    else:
        emp_str = '_noEMP'
        
    if p['use_Bulk_Formula']:
        qflux_str = '_bulk'
    else:
        qflux_str = ''
        
    if p['qnet_offset'] ==0:
        q_offset_str = ''
    else:
        q_offset_str = '_qoff%s'%p['qnet_offset']
        
    if p['iceMod']==1:
        qflux_str = '_iceModT'
    elif p['iceMod']==0:
        qflux_str = '_iceMod0'
        
    
    fnum = 9094#'0068' #9099
    p1 = 10#25 #10
    p2 = 25#28 #12
    nump = p2-p1
    met_data = 'NCEP_forcing_for_f%s_p%s-%s.nc' %(fnum, p1, p2)
    prof_data = 'float%s_%s_%s.nc' %(fnum, p1, p2)
    suffix = '%s_%sc%s%s%s%s_alpha%s%s' %(fnum, nump, ice_str, wind_str, emp_str, qflux_str, p['alpha'], q_offset_str)
    forcing, pwp_out = run_PWP(met_data=met_data, prof_data=prof_data, param_kwds=p, suffix=suffix, save_plots=True)
    
    return forcing, pwp_out  
    
    
def demo4(period='sum_win_2015'):
    
    """
    Example script of how to run the PWP model.
    This run is initialized with an late summer profile from float 9094 and forced by NCEP fluxes through till early winter. 
    
    test cases are for 2015 and 2016
    
    period == 'sum_win_2015', 'sum_win_2016', mix1, mix2
    
    """
    
    date_range_2015 = ['2015-Mar-03', '2015-Jul-15']
    date_range_2016 = ['2016-Feb-27', '2016-Jul-10']
    
    if period=='sum_win_2015':
        float_date_range = date_range_2015
        forcing_date_range = date_range_2015
    elif period=='sum_win_2016':
        float_date_range = date_range_2016
        forcing_date_range = date_range_2016
    elif period=='mix1':
        float_date_range = date_range_2015
        forcing_date_range = date_range_2016
    elif period=='mix2':
        float_date_range = date_range_2016
        forcing_date_range = date_range_2015
        
    
    
    p={}
    p['rkz']=1e-6
    p['dz'] = 1
    p['dt'] = 1.5
    p['max_depth'] = 500
    p['ice_ON'] = True
    p['winds_ON'] = True
    p['emp_ON'] = False
    p['alpha'] = 0.95
    p['dopt'] = 'pdens'
    p['fix_alpha'] = True
    p['mld_thresh'] = 0.01
    p['use_Bulk_Formula'] = True
    p['qnet_offset'] = -40 #W/m2

    
    if p['ice_ON']:
        ice_str = '' 
    else:
        ice_str = '_noICE'
        
    if p['winds_ON']:
        wind_str = ''
    else:
        wind_str = '_noWINDS'
        
    if p['emp_ON']:
        emp_str = ''
    else:
        emp_str = '_noEMP'
        
    if p['use_Bulk_Formula']:
        qflux_str = '_bulk'
    else:
        qflux_str = ''
        
    if p['qnet_offset'] ==0:
        q_offset_str = ''
    else:
        q_offset_str = '_qoff%s'%p['qnet_offset']
    
    fnum = 9094#'0068' #9099
    met_data = 'NCEP_forcing_for_f%s_%s-%s.nc' %(fnum, forcing_date_range[0], forcing_date_range[-1])
    prof_data = 'float%s_%s_%s.nc' %(fnum, float_date_range[0], float_date_range[-1])
    suffix = 'demo4_%s_float%s%s%s%s%s_alpha%s%s' %(period, fnum, ice_str, wind_str, emp_str, qflux_str, p['alpha'], q_offset_str)
    forcing, pwp_out = run_PWP(met_data=met_data, prof_data=prof_data, param_kwds=p, suffix=suffix, save_plots=True)
    
    return forcing, pwp_out, suffix 
    
    

def demo5(period='sum_win_2015'):
    
    """
    Like the demo4 but uses a summer profile from  float 9099
    
    """
    
    date_range_2015 = ['2015-Mar-02', '2015-Jul-14']
    date_range_2016 = ['2016-Feb-25', '2016-Jul-19']
    
    if period=='sum_win_2015':
        float_date_range = date_range_2015
        forcing_date_range = date_range_2015
    elif period=='sum_win_2016':
        float_date_range = date_range_2016
        forcing_date_range = date_range_2016
    elif period=='mix1':
        float_date_range = date_range_2015
        forcing_date_range = date_range_2016
    elif period=='mix2':
        float_date_range = date_range_2016
        forcing_date_range = date_range_2015
        
    
    
    p={}
    p['rkz']=1e-6
    p['dz'] = 1
    p['dt'] = 1.5
    p['max_depth'] = 500
    p['ice_ON'] = True
    p['winds_ON'] = True
    p['emp_ON'] = False
    p['alpha'] = 0.95
    p['dopt'] = 'pdens'
    p['fix_alpha'] = True
    p['mld_thresh'] = 0.01
    p['use_Bulk_Formula'] = True
    p['qnet_offset'] = -40 #W/m2

    
    if p['ice_ON']:
        ice_str = '' 
    else:
        ice_str = '_noICE'
        
    if p['winds_ON']:
        wind_str = ''
    else:
        wind_str = '_noWINDS'
        
    if p['emp_ON']:
        emp_str = ''
    else:
        emp_str = '_noEMP'
        
    if p['use_Bulk_Formula']:
        qflux_str = '_bulk'
    else:
        qflux_str = ''
        
    if p['qnet_offset'] ==0:
        q_offset_str = ''
    else:
        q_offset_str = '_qoff%s'%p['qnet_offset']
    
    fnum = 9099#'0068' #9099
    #met_data = 'NCEP_forcing_for_f%s_%s-%s.nc' %(fnum, forcing_date_range[0], forcing_date_range[-1])
    met_data = 'NCEP_forcing_for_f%s_%s-%s.nc' %(9094, '2015-Mar-03', '2015-Jul-15')
    prof_data = 'float%s_%s_%s.nc' %(fnum, float_date_range[0], float_date_range[-1])
    suffix = 'demo4_%s_float%s%s%s%s%s_alpha%s%s' %(period, fnum, ice_str, wind_str, emp_str, qflux_str, p['alpha'], q_offset_str)
    forcing, pwp_out = run_PWP(met_data=met_data, prof_data=prof_data, param_kwds=p, suffix=suffix, save_plots=True)
    
    return forcing, pwp_out, suffix 
    

def run_PWP(met_data, prof_data, param_kwds=None, overwrite=True, makeLivePlots=False, suffix='', save_plots=False):
    
    """
    This is the main controller function for the model. The flow of the algorithm
    is as follows:
        
        1) Set model parameters (see set_params function in PWP_helper.py).
        2) Read in forcing and initial profile data.
        3) Prepare forcing and profile data for model run (see prep_data in PWP_helper.py).
            3.1) Interpolate forcing data to prescribed time increments.
            3.2) Interpolate profile data to prescribed depth increments.
            3.3) Initialize model output variables.
        4) Iterate the PWP model specified time interval:
            4.1) apply winds, heat and salt fluxes
            4.2) rotate, adjust to wind, rotate
            4.3) apply bulk Richardson number mixing
            4.4) apply gradient Richardson number mixing
            4.5) apply drag associated with internal wave dissipation
            4.5) apply diapycnal diffusion
        5) Save results to output file
    
    Input:
    met_data -  path to netCDF file containing forcing/meterological data. This file must be in the
                input_data/ directory.
                
                The data fields should include 'time', 'sw', 'lw', 'qlat', 'qsens', 'tx',
                'ty', and 'precip'. These fields should store 1-D time series of the same
                length.
                
                The model expects positive heat flux values to represent ocean warming. The time
                data field should contain a 1-D array representing fraction of day. For example,
                for 6 hourly data, met_data['time'] should contain a number series that increases
                in steps of 0.25, such as np.array([1.0, 1.25, 1.75, 2.0, 2.25...]).
                
                See https://github.com/earlew/pwp_python#input-data for more info about the
                expected intput data.
    
    prof_data - path to netCDF file containing initial profile data. This must be in input_data/ directory.
                Fields should include: ['z', 't', 's', 'lat']. These represent 1-D vertical profiles of temperature,
                salinity and density. lat can (should?) be a float.
                
                See https://github.com/earlew/pwp_python#input-data for more info about the
                expected intput data
    
    overwrite - controls the naming of output file. If True, the same filename is used for
                every model run. If False, a unique time_stamp is generated and appended
                to the file name. Default is True.
    
    makeLivePlots - if True, the code will generate live plots of mixed layer properties at
                each time step (makes the code run a lot SLOWER). Default is False
    
    suffix - string to add to the end of filenames. e.g. suffix = 'nodiff' leads to 'pwp_out_nodiff.nc.
            default is an empty string ''.
    
    save_plots -this gets passed on to the makeSomePlots() function in the PWP_helper. If True, the code
                saves the generated plots. Default is False.
    
    param_kwds -dict containing keyword arguments for set_params function. See PWP_helper.set_params()
                for more details. If None, default parameters are used. Default is None.
    
    Output:
    
    forcing, pwp_out = PWP.run()
    
    forcing: a dictionary containing the interpolated surface forcing.
    pwp_out: a dictionary containing the solutions generated by the model.
    
    This script also saves the following to file:
    
    'pwp_output.nc'- a netCDF containing the output generated by the model.
    'pwp_output.p' - a pickle file containing the output generated by the model.
    'forcing.p' - a pickle file containing the (interpolated) forcing used for the model run.
    If overwrite is set to False, a timestamp will be added to these file names.
    
    ------------------------------------------------------------------------------
    There are two ways to run the model:
    1.  You can run the model by typing "python PWP.py" from the bash command line. This
        will initiate this function with the set defaults. Typing "%run PWP" from the ipython
        command line will do the same thing.
    
    2.  You can also import the module then call the run() function specifically. For example,
        >> import PWP
        >> forcing, pwp_out = PWP.run()
        Alternatively, if you want to change the defaults...
        >> forcing, pwp_out = PWP.run(met_data='new_forcing.nc', overwrite=False, diagnostics=False)
    
    This is a more interactive approach as it provides direct access to all of the model's
    subfunctions.
    
    """
    
    #close all figures
    plt.close('all')
    
    #start timer
    t0 = timeit.default_timer()
    
    ## Get surface forcing and profile data
    met_dset = xr.open_dataset('input_data/%s' %met_data)
    prof_dset = xr.open_dataset('input_data/%s' %prof_data)
    
    if 'start_date' in prof_dset.attrs and 'end_date' in prof_dset.attrs:
        
        print('=============================================')
        print("Initializing with float %s..." %prof_dset.float_num)
        print("Float starting location: %.2fE, %.2fN" %(prof_dset['lon'][0].values, prof_dset['lat'][0].values))
        print("Start time: %s" %(prof_dset.start_date))
        print("Approximate end time: %s" %(prof_dset.end_date))
        print('=============================================')
    
    ## get model parameters and constants (read docs for set_params function)
    try:
        lat = prof_dset['lat'].values[0] #needed to compute internal wave dissipation
    except IndexError:
        lat = prof_dset['lat'].values
    
    if param_kwds is None:
        params = set_params(lat=lat)
    else:
        param_kwds['lat'] = lat
        params = set_params(**param_kwds)
    
    ## prep forcing and initial profile data for model run (see prep_data function for more details)
    forcing, pwp_out, params = prep_data(met_dset, prof_dset, params)
    
    prof_dset.close()
    met_dset.close()
    
    # plot forcing
    makeSomePlots(forcing, pwp_out, save_plots=True, justForcing=True)
    
    ## run the model
    pwp_out = PWP.pwpgo(forcing, params, pwp_out, makeLivePlots)
    
    #check timer
    tnow = timeit.default_timer()
    t_elapsed  = (tnow - t0)
    print("Time elapsed: %i minutes and %i seconds" %(np.floor(t_elapsed/60), t_elapsed%60))
    
    
    ## save output as pickle files
    if overwrite:
        time_stamp = ''
    else:
        #use unique time stamp
        time_stamp = datetime.now().strftime("_%Y%m%d_%H%M")

    if len(suffix)>0 and suffix[0] != '_':
        suffix = '_%s' %suffix
     
    output_fpath = "output/pwp_output%s%s.p" %(suffix, time_stamp)
    forcing_fpath = "output/forcing%s%s.p" %(suffix, time_stamp)
    
    pickle.dump(pwp_out, open(output_fpath, 'wb'))
    pickle.dump(forcing, open(forcing_fpath, 'wb'))
    

    # # save output as netCDF file
    # output_fpath = "output/pwp_output%s%s.nc" %(suffix, time_stamp)
    # forcing_fpath = "output/forcing%s%s.nc" %(suffix, time_stamp)
    #
    # pwp_out2 = phf.save2nc(pwp_out, output_fpath, dt_save=params['dt_save'])
    # forcing2 = phf.save2nc(forcing, forcing_fpath, dt_save=params['dt_save'], type='forc')
    
    #debug_here()
    ## do analysis of the results
    makeSomePlots(forcing, pwp_out, zlim=params['plot_zlim'], suffix=suffix, save_plots=save_plots)
    
    return forcing, pwp_out   


def set_params(lat, dt=3., dz=1., max_depth=100., mld_thresh=1e-4, dt_save=1, alpha=-999, h_i0=0.0, rkz=0., diff_zlim=5000, plot_zlim=500, qnet_offset=0., dopt='dens0', use_Bulk_Formula=False, fix_alpha=False, ice_ON=False, winds_ON=True, emp_ON=True, drag_ON=True, iceMod=1):
                
    
    """
                
    TODO: combine all these input into dict
    
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
    
    
    rkz: background vertical diffusion (m**2/s). [0.]
    
    
    OUTPUT is dict with fields containing the above variables plus the following:
    dt_d: time increment (dt) in units of days
    g: acceleration due to gravity [9.8 m/s^2]
    cpw: specific heat of water [4183.3 J/kgC]
    f: coriolis term (rad/s). [sw.f(lat)]
    ucon: coefficient of inertial-internal wave dissipation (s^-1) [0.1*np.abs(f)]
                
    """
    
    rb=0.65 #rb: critical bulk richardson number. [0.65]
    rg=0.25 #rg: critical gradient richardson number. [0.25]
    beta1=0.6 #beta1: longwave extinction coefficient (meters). [0.6]
    beta2=20.0 #beta2: shortwave extinction coefficient (meters). [20] 
    
    params = {}
    params['dt'] = 3600.0*dt
    params['dt_d'] = params['dt']/86400.
    params['dz'] = dz
    params['dt_save'] = int(dt_save)
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
    params['diff_zlim'] = diff_zlim
    params['plot_zlim'] = plot_zlim
    
    params['ice_ON'] = ice_ON
    params['winds_ON'] = winds_ON
    params['emp_ON'] = emp_ON
    params['drag_ON'] = drag_ON
    
    
    params['iceMod'] = iceMod #
    params['alpha'] = alpha #sea ice concentration 
    params['fix_alpha'] = fix_alpha #if fix_alpha=True, use params['alpha'] instead of ice conc. from forcing
    params['h_i0'] = h_i0 #initial ice thickness
    
    params['qnet_offset'] = qnet_offset #arbitrary offset to the net atmospheric heat flux.
    
    params['dens_option'] = dopt # 'dens', 'dens0' or 'pdens' 
    """
    NOTES about density option:
    default is dens0 because that was the original specification.
    - 'dens0' calls dens0() from the seawater module which computes density assuming s(p=0) and t(p=0). Note that this is not potential density.
    - 'pdens' calls pden() and computes potential density, i.e. density using potential temp with a surface reference
    - 'dens' cals dens(), which uses the full density equation.

    'dens0' is the simpliest/fastest option but can produce weak density inversions in weakly statified but otherwise stable water columns.
    
    """
    
    
    params['use_Bulk_Formula'] = use_Bulk_Formula #compute surface fluxes using bulk formulae (True) or simply use what's provided (False)
    
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
            ['z', 't', 's', 'lat'], where 't' and 's' represent 1-D vertical profiles of temperature,
            salinity. 
    
            NOTE: Code now accepts two column arrays. Second column will be treated as observed
            profile at the end of run.
    
            update: code now accepts a 'ps' vector which represents a 1-D passive scalar
            
    params: dictionary-like object with fields defined by set_params function
    
    OUTPUT:
    
    forcing: dictionary with interpolated surface forcing data. 
    pwp_out: dictionary with initialized variables to collect model output.
    """
    
    import warnings
    
    #create new time vector with time step dt_d
    #time_vec = np.arange(met_dset['time'][0], met_dset['time'][-1]+params['dt_d'], params['dt_d']) 
    time_vec = np.arange(met_dset['time'][0], met_dset['time'][-1], params['dt_d']) 
    tlen = len(time_vec)
    
    #interpolate surface forcing data to new time vector
    from scipy.interpolate import interp1d
    forcing = {} 
    for vname in met_dset:
        p_intp = interp1d(met_dset['time'], met_dset[vname], axis=0)
        forcing[vname] = p_intp(time_vec)
        
    if 'dtime' in met_dset:
        from datetime import timedelta
        dt_model = timedelta(seconds=params['dt'])
        dtime2 = np.arange(met_dset['dtime'][0].values, met_dset['dtime'][-1].values, dt_model)
        forcing['dtime'] = dtime2
        
    #adjust skin temperature if it exists
    if 'skt' in list(forcing.keys()):
        #if met_dset['skt'].attrs['units'] == 'degK':
        forcing['skt'] = forcing['skt'] - 273.15
        #forcing['skt'] = 0.75*forcing['skt']
        # forcing['skt'].attrs['units'] = 'degC'
    else:
        forcing['skt'] = np.zeros(len(forcing['sw']))*np.nan
        
        
    #if icec doesn't exit, add it as dummy variable
    if 'icec' not in list(forcing.keys()):
        forcing['icec'] = np.zeros(len(forcing['sw']))*np.nan
        
        
    #arbitrarily adjust (tune) fluxes
    # forcing['qsens'][:] = 0.0
    # forcing['sw'] = 0.75*forcing['sw']
    # # forcing['lw'] = 1.2*forcing['lw']
    # # forcing['qlat'] = 1.2*forcing['qlat']
    # forcing['precip'] = 0.5*forcing['precip']
    # print "WARNING: fluxes were adjusted!!!"
            
        
    #interpolate E-P to dt resolution (not sure why this has to be done separately)
    evap_intp = interp1d(met_dset['time'], met_dset['qlat'], axis=0, kind='nearest', bounds_error=False)
    evap = (0.03456/(86400*1000))*evap_intp(np.floor(time_vec)) #(meters per second)
    emp = evap-forcing['precip']
    emp[np.isnan(emp)] = 0.
    forcing['emp'] = emp  
    
    if params['emp_ON'] == False:
        print("E-P is turned OFF.")
        forcing['emp'][:] = 0.0
          
    
    #define F_in and F_out 
    #for sw, lw, qlat and qsens and positive values should mean ocean warming
    #for the purpose of the PWP code, we flip the sign of F_out
    
    forcing['F_in'] = forcing['sw'] #heat flux into ocean
    forcing['F_out'] = -(forcing['lw'] + forcing['qlat'] + forcing['qsens']) 
    
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
    
    #compute absorption and incoming radiation (function defined in PWP.py)
    absrb = PWP.absorb(params['beta1'], params['beta2'], zlen, params['dz']) #(units unclear)
    
    #check for numeric stability. This relates to the diffusion equation
    dstab = params['dt']*params['rkz']/params['dz']**2 #courant number  
    if dstab > 0.5:
        print("WARNING: unstable CFL condition for diffusion! dt*rkz/dz**2 > 0.5.")
        print("To fix this, try to reduce the time step or increase the depth increment.")
        inpt = eval(input("Proceed with simulation? Enter 'y'or 'n'. "))
        if inpt is 'n':
            raise ValueError("Please restart PWP.m with a larger dz and/or smaller dt. Exiting...")
        
    forcing['absrb'] = absrb
    params['dstab'] = dstab
    
    #check depth resolution of profile data
    prof_z = prof_dset['z']; max_z = params['max_depth']
    prof_incr = np.diff(prof_dset['z'][prof_z<=max_z]).mean()
    if params['dz'] < prof_incr/5.:
        message = "Specified depth increment (%s m), is much smaller than mean profile resolution (%s m)." %(params['dz'], prof_incr)
        #warnings.warn(message)
        print(message)
        
        
        # inpt = input("Depth increment, dz, is much smaller than profile resolution. Is this okay? (Enter 'y'or 'n')")
        # if inpt is 'n':
        #     raise ValueError("Please restart PWP.m with a new dz >= %s Exiting..." %prof_incr/5.)
    
    #debug_here()
    #interpolate profile data to new z-coordinate
    from scipy.interpolate import InterpolatedUnivariateSpline 
    for vname in prof_dset:
        if vname == 'lat' or vname=='lon' or vname=='p' or vname=='z':
            continue
        else:
            
            vble = prof_dset[vname].values
            
            if vble.ndim==1:
                vble = vble[:, np.newaxis]
            
            assert vble.ndim>1, "profile variable should be 2-D"
            
            #first strip nans
            not_nan = np.logical_not(np.isnan(vble[:,0]))
            if vname=='ps' and np.all(np.isnan(vble[:,0])):
                init_prof[vname] = np.zeros(init_prof['z'].shape)
            
            else:
                
                #indices = np.arange(len(vble[:,0]))
                #p_intp = interp1d(prof_dset['z'], prof_dset[vname], axis=0, kind='linear', bounds_error=False)
                #interp1d doesn't work here because it doesn't extrapolate. Can't have Nans in interpolated profile
                p_intp = InterpolatedUnivariateSpline(prof_dset['z'][not_nan], vble[not_nan, 0], k=1)
                init_prof[vname] = p_intp(init_prof['z'])    
        
    #get profile variables
    temp0 = init_prof['t'] #initial profile temperature
    sal0 = init_prof['s'] #intial profile salinity
    
    #get passive scalar if any
    if 'ps' in prof_dset:
        ps0 = init_prof['ps'] #initial profile temperature of passive scalar
    else:
        ps0 = np.zeros(temp0.shape)
    
    print("Using %s density option." %params['dens_option'])
    
    dens0 = PWP.getDensity(sal0, temp0, init_prof['z'], dopt=params['dens_option'])
        
    if any(np.diff(dens0)<0):
        message = "Warning!!! Initial density profile has instabilities..."
        # print(np.diff(dens0))
        #warnings.warn(message)
        print(message)
        
    sal0, temp0, ps0, dens0  = PWP.local_stir(init_prof['z'], sal0, temp0, ps0, params['dens_option'])
    
    #initialize variables for output
    #Todo: set time resolution of output file
    pwp_out = {}
    pwp_out['time'] = time_vec
    pwp_out['dt'] = params['dt']
    pwp_out['dz'] = params['dz']
    pwp_out['lat'] = params['lat']
    pwp_out['z'] = init_prof['z']
    
    #tlen = int(np.floor(tlen/params['dt_save']))
    arr_sz = (zlen, tlen)
    pwp_out['temp'] = np.zeros(arr_sz)*np.nan
    pwp_out['sal'] = np.zeros(arr_sz)*np.nan
    pwp_out['dens'] = np.zeros(arr_sz)*np.nan
    pwp_out['uvel'] = np.zeros(arr_sz)
    pwp_out['vvel'] = np.zeros(arr_sz)
    pwp_out['ps'] = np.zeros(arr_sz)*np.nan
    pwp_out['mld'] = np.zeros((tlen,))*np.nan
    
    #use temp, sal, dens and passive tracer profile data for the first time step
    pwp_out['sal'][:,0] = sal0
    pwp_out['temp'][:,0] = temp0
    pwp_out['dens'][:,0] = dens0
    pwp_out['ps'][:,0] = ps0
    
    #find initial ml index
    mld_idx = np.flatnonzero(dens0-dens0[0]>params['mld_thresh'])[0]
    pwp_out['mld'][0] = pwp_out['z'][mld_idx] 
    
    #if final observed profile is available, save it
    if prof_dset['t'].ndim == 2:
        pwp_out['obs_zlvl'] = prof_dset['z'].values
        pwp_out['sal_f'] = prof_dset['s'][:,1].values
        pwp_out['temp_f'] = prof_dset['t'][:,1].values
        

    #create variables for sea ice
    pwp_out['surf_ice_temp'] = np.nan*np.zeros((tlen,))
    pwp_out['ice_thickness'] = np.zeros((tlen,))
    
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
    ps = pwp_out['ps']
    z = pwp_out['z']


    #plot depth int. KE and momentum
    # plt.figure(num=101)
    #
    # plt.subplot(211)
    # plt.plot(time[n]-time[0], np.trapz(0.5*dens[:,n]*(uvel[:,n]**2+vvel[:,n]**2)), 'b.')
    # plt.grid(True)
    # if n==1:
    #     plt.title('Depth integrated KE')
    #
    # plt.subplot(212)
    # plt.plot(time[n]-time[0], np.trapz(dens[:,n]*np.sqrt(uvel[:,n]**2+vvel[:,n]**2)), 'b.')
    # plt.grid(True)
    # plt.pause(0.05)
    # plt.subplots_adjust(hspace=0.35)

    # #debug_here()
    # if n==1:
    #     plt.title('Depth integrated Mom.')
    #     #plt.get_current_fig_manager().window.wm_geometry("400x600+20+40")
    
    #plot T,S and U,V
    plt.figure(num=102, figsize=(12.5,7))
    ax1 = plt.subplot2grid((1,5), (0, 0), colspan=2)
    ax1.plot(uvel[:,n], z, 'b', label='uvel')
    ax1.plot(vvel[:,n], z, 'r', label='vvel')
    ax1.invert_yaxis()
    # xlims = ax1.get_xlim()
    # xticks = np.round(np.linspace(xlims[0], xlims[1], 4), 1)
    # ax1.set_xticks(xticks)
    ax1.set_xlim(-0.06, 0.06)
    ax1.grid(True)
    ax1.legend(loc=3)    

    ax2 = plt.subplot2grid((1,5), (0, 2), colspan=1)
    ax2.plot(temp[:,n], z, 'b')
    ax2.grid(True)
    ax2.set_xlabel('Temp.')
    ax2.invert_yaxis()
    xlims = ax2.get_xlim()
    xticks = np.round(np.linspace(xlims[0], xlims[1], 3), 1)
    ax2.set_xticks(xticks)

    ax3 = plt.subplot2grid((1,5), (0, 3), colspan=1)
    ax3.plot(sal[:,n], z, 'b')
    ax3.set_xlabel('Salinity')
    ax3.grid(True)
    ax3.invert_yaxis()
    xlims = ax3.get_xlim()
    xticks = np.round(np.linspace(xlims[0], xlims[1], 3), 1)
    ax3.set_xticks(xticks)
    
    ax4 = plt.subplot2grid((1,5), (0, 4), colspan=1)
    ax4.plot(ps[:,n], z, 'b')
    ax4.set_xlabel('Passive Scalar')
    ax4.grid(True)
    ax4.invert_yaxis()
    xlims = ax4.get_xlim()
    xticks = np.round(np.linspace(xlims[0], xlims[1], 3), 0)
    ax4.set_xticks(xticks)

    plt.subplots_adjust(wspace=0.4)
    plt.pause(0.05)

    plt.show()
    
    
def formatDates(ax):
    
    import matplotlib.dates as mdates
    # years_md = mdates.YearLocator()   # every year
    months_md = mdates.MonthLocator()  # every month
    days_md15 = mdates.DayLocator([1,15])  # every 15 days
    days_md1 = mdates.DayLocator()  # every 1 day
    dateFmt = mdates.DateFormatter('%b-%d-%Y')
    #yearsFmt = mdates.DateFormatter('%Y')
    
    ax.xaxis.set_major_locator(days_md15)
    ax.xaxis.set_major_formatter(dateFmt)
    ax.xaxis.set_minor_locator(days_md1)   
    
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)


def makeSomePlots(forcing, pwp_out, zlim=500, save_plots=False, suffix='', justForcing=False, showPlots=True):
    
    """
    TODO: add doc file
    Function to make plots of the results once the model iterations are complete.
    
    """
    
    if len(suffix)>0 and suffix[0] != '_':
            suffix = '_%s' %suffix
    
    ## Plot surface fluxes
    fig0, axes = plt.subplots(3,1, sharex=True, figsize=(7.5,9))
    
    if 'dtime' in forcing:
        tvec = forcing['dtime']
    else:
        tvec = pwp_out['time']
    
    # if time_vec is None:
    #     tvec = pwp_out['time']
    # else:
    #     tvec = time_vec
    #
    
    axes = axes.flatten()
    ## plot surface heat flux
    axes[0].plot(tvec, forcing['lw'], label='$Q_{lw}$')
    axes[0].plot(tvec, forcing['qlat'], label='$Q_{lat}$')
    axes[0].plot(tvec, forcing['qsens'], label='$Q_{sens}$')
    axes[0].plot(tvec, forcing['sw'], label='$Q_{sw}$')
    axes[0].hlines(0, tvec[0], tvec[-1], linestyle='-', color='0.3')
    axes[0].plot(tvec, forcing['F_in']-forcing['F_out'], ls='-', lw=2, color='k', label='$Q_{net}$')   
    axes[0].set_ylabel('Heat flux (W/m2)')
    axes[0].set_title('Heat flux into ocean')
    axes[0].grid(True)
    axes[0].set_xlim()
    
    ##plot wind stress
    axes[1].plot(tvec, forcing['tx'], label=r'$\tau_x$')
    axes[1].plot(tvec, forcing['ty'], label=r'$\tau_y$')
    axes[1].hlines(0, tvec[0], tvec[-1], linestyle='--', color='0.3')
    axes[1].set_ylabel('Wind stress (N/m2)')
    axes[1].set_title('Wind stress')
    axes[1].grid(True)
    axes[1].legend(loc=0, fontsize='medium')
    
    
    ## plot freshwater forcing
    emp_mmpd = forcing['emp']*1000*3600*24 #convert to mm per day
    axes[2].plot(tvec, emp_mmpd, label='E-P')
    axes[2].hlines(0, tvec[0], tvec[-1], linestyle='--', color='0.3')
    axes[2].set_ylabel('Freshwater forcing (mm/day)')
    axes[2].set_title('Freshwater forcing')
    axes[2].grid(True)
    axes[2].legend(loc=0, fontsize='medium')
    axes[2].set_xlabel('Time (days)')
    
    if 'dtime' in forcing:
        formatDates(axes[2])
    
    ## plot ice conc. and surf temp.
    fig1, axes = plt.subplots(2,1, sharex=True, figsize=(7.5,9))
    axes = axes.flatten()
    axes[0].plot(tvec, forcing['icec'], label='Forced')
    axes[0].set_ylabel('Ice cover percentage (%)')
    axes[0].set_title('Ice concentration')
    axes[0].grid(True)
    
    axes[1].plot(tvec, forcing['skt'])
    axes[1].set_ylabel('Temperature (C)')
    axes[1].set_title('Skin temperature')
    axes[1].grid(True)
    axes[1].set_xlabel('Time (days)')

    if 'dtime' in forcing:
        formatDates(axes[1])
        
        
    
    
    # plot surface variables used to compute bulk fluxes
    fig, axes = plt.subplots(3,1, sharex=True, figsize=(7.5,9))
    axes = axes.flatten()
    axes[0].plot(tvec, forcing['atemp2m']-273.15, label='2m Temp')
    axes[0].plot(tvec, forcing['skt'], label='Skin Temp')
    axes[0].legend(loc=0, fontsize='medium')
    axes[0].set_ylabel('Temperature ($^{\circ}$C)')
    axes[0].set_title('Near-surface temperature')
    axes[0].grid(True)
    
    axes[1].plot(tvec, forcing['u10m'], label='u-wind')
    axes[1].plot(tvec, forcing['v10m'], label='v-wind')
    axes[1].legend(loc=0, fontsize='medium')
    axes[1].set_ylabel('Wind speed')
    axes[1].set_title('10 meter winds (m/s)')
    axes[1].grid(True)
    
    axes[2].plot(tvec, forcing['shum2m']*1000)
    axes[2].set_ylabel('Specific Humidity (g/kg)')
    axes[2].set_title('2m specific humidity')
    axes[2].grid(True)
    
    if 'dtime' in forcing:
        formatDates(axes[2])
        
    if save_plots:     
        fig.savefig('plots/surface_state%s.pdf' %suffix, bbox_inches='tight')

    if justForcing:
        plt.show()
        plt.pause(1)
        return
        
    #add true ice conc. to above plot
    axes[0].plot(tvec, pwp_out['alpha_true'], label='Actual')
    axes[0].legend(loc=0, fontsize='medium')
    
    
    
    if save_plots:     
        fig0.savefig('plots/surface_forcing%s.pdf' %suffix, bbox_inches='tight')
        fig1.savefig('plots/ice_conc_temp%s.pdf' %suffix, bbox_inches='tight')
        
        plt.close('all')
    
    
    ## Plot computed heat fluxes
    fig2, axes = plt.subplots(3,1, sharex=True, figsize=(7.5,9))
    axes = axes.flatten()
    ## plot atmosphere ocean flux
    axes[0].plot(tvec, pwp_out['F_lw_ao'], label='$Q_{lw}^{ocn}$')
    axes[0].plot(tvec, pwp_out['F_lat_ao'], label='$Q_{lat}^{ocn}$')
    axes[0].plot(tvec, pwp_out['F_sens_ao'], label='$Q_{sens}^{ocn}$')
    axes[0].plot(tvec, (1-pwp_out['alpha_true'])*forcing['sw'], label='$Q_{sw}^{ocn}$')
    axes[0].hlines(0, tvec[0], tvec[-1], linestyle='-', color='0.3')
    axes[0].plot(tvec, pwp_out['F_net_ao'], ls='-', lw=2, color='k', label='$Q_{net}^{ocn}$')
    axes[0].set_ylabel('Heat flux (W/m2)')
    axes[0].set_title('Computed Atmosphere-ocean heat flux')
    axes[0].grid(True)
    #axes[0].set_ylim(-500,300)

    axes[0].legend(loc=0, ncol=2, fontsize='medium')

    ## plot atmosphere ice flux
    axes[1].plot(tvec, pwp_out['F_lw_ai'], label='$Q_{lw}^{ice}$')
    axes[1].plot(tvec, pwp_out['F_lat_ai'], label='$Q_{lat}^{ice}$')
    axes[1].plot(tvec, pwp_out['F_sens_ai'], label='$Q_{sens}^{ice}$')
    axes[1].plot(tvec, pwp_out['alpha_true']*forcing['sw'], label='$Q_{sw}^{ice}$')
    axes[1].hlines(0, tvec[0], tvec[-1], linestyle='-', color='0.3')
    axes[1].plot(tvec, pwp_out['F_net_ai'], ls='-', lw=2, color='k', label='$Q_{net}^{ice}$')
    axes[1].set_ylabel('Heat flux (W/m2)')
    axes[1].set_title('Computed atmosphere-ice heat flux')
    axes[1].grid(True)
    #axes[0].set_ylim(-500,300)

    axes[1].legend(loc=0, ncol=2, fontsize='medium')
    
    ## plot atmosphere ice flux
    axes[2].plot(tvec, pwp_out['F_lw_ai']+pwp_out['F_lw_ao'], label='$Q_{lw}$')
    axes[2].plot(tvec, pwp_out['F_lat_ai']+pwp_out['F_lat_ao'], label='$Q_{lat}$')
    axes[2].plot(tvec, pwp_out['F_sens_ai']+pwp_out['F_sens_ao'], label='$Q_{sens}$')
    axes[2].plot(tvec, forcing['sw'], label='$Q_{sw}$')
    axes[2].hlines(0, tvec[0], tvec[-1], linestyle='-', color='0.3')
    axes[2].plot(tvec, pwp_out['F_net_ai'] + pwp_out['F_net_ai'], ls='-', lw=2, color='k', label='$Q_{net}$')
    axes[2].set_ylabel('Heat flux (W/m2)')
    axes[2].set_title('Total surface heat flux')
    axes[2].grid(True)
    #axes[0].set_ylim(-500,300)

    axes[2].legend(loc=0, ncol=2, fontsize='medium')
    
    if 'dtime' in forcing:
        formatDates(axes[2])
        
    
    if save_plots:     
        fig2.savefig('plots/computed_surface_forcing%s.pdf' %suffix, bbox_inches='tight')
        
        plt.close('all')
    
    
    
    
    #plot summary of ML evolution
    
    #get a smoothed MLD time series
    from scipy.signal import savgol_filter
    dt = pwp_out['time'][1]-pwp_out['time'][0]
    mld_exact2 = pwp_out['mld_exact2'].copy()
    mld_exact2_sm = np.zeros(len(mld_exact2))*np.nan
    nan_i = np.isnan(mld_exact2)
    wlen = (5*(1./dt)+1)
    if wlen%2==0: wlen+=1
    mld_exact2_sm[~nan_i] = savgol_filter(mld_exact2[~nan_i], window_length=wlen, polyorder=1, mode='nearest')
    
    ##plot temp and sal change over time
    fig, axes = plt.subplots(3,1, sharex=True, figsize=(8,8))
    vble = ['temp', 'sal', 'ps']
    units = ['$^{\circ}$C', 'PSU', '$\mu$mol/kg']
    #cmap = custom_div_cmap(numcolors=17)
    clim = ([-1.5, 2.0], [33.75, 34.75], [200, 325])
    cmap = [plt.cm.coolwarm, plt.cm.RdYlGn_r, plt.cm.rainbow]
    for i in range(3):
        ax = axes[i]
        clvls = np.linspace(clim[i][0], clim[i][1], 21)
        im = ax.contourf(pwp_out['time'], pwp_out['z'], pwp_out[vble[i]], clvls, cmap=cmap[i], extend='both')
        ax.plot(tvec[1:], mld_exact2_sm[1:], '-k')
        ax.set_ylabel('Depth (m)')
        ax.set_ylim(0,zlim)
        ax.set_title('Evolution of ocean %s (%s)' %(vble[i], units[i]))
        ax.invert_yaxis()   
        cb = plt.colorbar(im, ax=ax, format='%.2f')
     
    ax.set_xlabel('Days') 
    
    if 'dtime' in forcing:
        formatDates(ax)
    
    if save_plots:     
        plt.savefig('plots/sal_temp_ztseries_%s.png' %suffix, bbox_inches='tight') 
        plt.close()
    
    
    ## plot initial and final T-S profiles (TODO: add actual T,S profiles)
    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    
    plt.figure()
    host = host_subplot(111, axes_class=AA.Axes)
    par1 = host.twiny() #par for parasite axis
    host.set_ylabel("Depth (m)")
    host.set_xlabel("Temperature ($^{\circ}$C)")
    par1.set_xlabel("Salinity (PSU)")
    
    host.set_ylim(0, max(pwp_out['z']))
    host.invert_yaxis()
    
    p1, = host.plot(pwp_out['temp'][:,0], pwp_out['z'], '--r', label='$T_i$')
    host.plot(pwp_out['temp'][:,-1], pwp_out['z'], '-r', label='$T_f$')
    p2, = par1.plot(pwp_out['sal'][:,0], pwp_out['z'], '--b', label='$S_i$')
    par1.plot(pwp_out['sal'][:,-1], pwp_out['z'], '-b', label='$S_f$')
    host.grid(True)
    
    host.legend(loc=3, ncol=2)
    #par1.legend(loc=3)
    
    #if observed final profiles are available, plot them:
    # if 'sal_f' in pwp_out.keys() and 'temp_f' in pwp_out.keys():
    #     host.plot(pwp_out['temp_f'], pwp_out['obs_zlvl'], '-o', color='r', label='$T_{obs}$')
    #     par1.plot(pwp_out['sal_f'], pwp_out['obs_zlvl'], '-o', color='b', label='$S_{obs}$')
    #
    
    
    host.axis["bottom"].label.set_color(p1.get_color())
    host.axis["bottom"].major_ticklabels.set_color(p1.get_color())
    host.axis["bottom"].major_ticks.set_color(p1.get_color())

    par1.axis["top"].label.set_color(p2.get_color())
    par1.axis["top"].major_ticklabels.set_color(p2.get_color())
    par1.axis["top"].major_ticks.set_color(p2.get_color())

    
    if save_plots:     
        plt.savefig('plots/initial_final_TS_profiles%s.pdf' %suffix, bbox_inches='tight')
        plt.close()    
        
    ## plot ice growth and ice temp
    fig, axes = plt.subplots(2,1, figsize=(7.5,9))
    axes[0].plot(tvec, pwp_out['ice_thickness'], '-')
    axes[0].set_ylabel('Ice thickness (m)')
    #axes[0].xlabel('Time (days)')
    axes[0].set_title('Ice thickness')
    axes[0].grid(True)
    
    if 'dtime' in forcing:
        formatDates(axes[0])
    
    axes[1].plot(tvec, pwp_out['surf_ice_temp'], '-')
    axes[1].set_ylabel('Temperature (C)')
    axes[1].set_xlabel('Time (days)')
    axes[1].set_title('Ice surface temperature')
    axes[1].grid(True)
    
    if 'dtime' in forcing:
        formatDates(axes[1])
    
    if save_plots:     
        plt.savefig('plots/ice_temp_thickness%s.pdf' %suffix, bbox_inches='tight')
        plt.close()
     
    plt.figure(figsize=(6,6.5))
    ax = plt.gca()
    im = ax.scatter(pwp_out['ice_thickness'], pwp_out['surf_ice_temp'], s=10, c=tvec, cmap=plt.cm.rainbow, lw=0)
    ax.plot(pwp_out['ice_thickness'], pwp_out['surf_ice_temp'], '-', lw=0.5, c='0.5')
    ax.set_ylabel('Temperature (C)')
    ax.set_xlabel('Ice thickness (m)')
    ax.set_title('Ice surface temperature')
    cb = plt.colorbar(im)
    cb.set_label('Time (days)')
    ax.grid(True)
    
    if save_plots:     
        plt.savefig('plots/ice_temp_thickness_v2%s.pdf' %suffix, bbox_inches='tight')
        plt.close()
    
        
    ## plot OCEAN-ICE heat flux and ATM-ICE heat flux
    F_oi = np.ma.masked_invalid(pwp_out['F_oi'])
    F_oi_arr = np.ma.filled(F_oi, 0)
    F_ocean_net = -F_oi_arr+np.ma.masked_invalid(pwp_out['F_ao'])
    
    
    F_ai_sm = savgol_filter(pwp_out['F_ai'], window_length=wlen, polyorder=1, mode='nearest')
    F_i_sm = savgol_filter(pwp_out['F_i'], window_length=wlen, polyorder=1, mode='nearest')
    F_oi_sm = savgol_filter(F_oi_arr, window_length=wlen, polyorder=1, mode='nearest')
    F_ocean_net_sm = savgol_filter(F_ocean_net, window_length=wlen, polyorder=1, mode='nearest')
    
    
    plt.figure()
    plt.subplot(111)
    plt.plot(tvec, pwp_out['F_ai'], lw=0.5, alpha=0.25, c='b')
    plt.plot(tvec, F_ai_sm, '-', label='$F_{ai}$', lw=2, c='b')
    
    plt.plot(tvec, pwp_out['F_i'], lw=0.5, alpha=0.25, c='g')
    plt.plot(tvec, F_i_sm, label='$F_i$', lw=2, c='g')
    
    plt.plot(tvec, pwp_out['F_oi'], lw=0.5, alpha=0.5, c='tomato')
    plt.plot(tvec, F_oi_sm, label='$F_{oi}$', lw=2, c='tomato')

    plt.plot(tvec,  F_ocean_net, lw=0.5, alpha=0.25, c='magenta')
    plt.plot(tvec,  F_ocean_net_sm, label='$F_{ao}$ - $F_{oi}$', lw=2, c='magenta')
    
    plt.hlines(0, tvec[0], tvec[-1])
    plt.ylim(-400, 200)
    #plt.ylim(-1.5*np.abs(pwp_out['F_atm'].max()), 1.5*np.abs(pwp_out['F_atm'].max()))
    plt.xlabel('Time (days)')
    plt.ylabel('Heat flux (W/m2)')
    plt.legend(loc=0, fontsize='small')
    plt.grid(True)
    
    if 'dtime' in forcing:
        formatDates(plt.gca())
    
    print("mean ocean warming: %.5f W/m2" %F_ocean_net.mean())
    
    if save_plots:     
        plt.savefig('plots/ice_ocean_fluxes%s.pdf' %suffix, bbox_inches='tight')
        plt.close()
        
    #plot change in MLD  
    plt.figure()
    plt.subplot(111)
    plt.plot(tvec[1:], pwp_out['mld'][1:], label='model approx.')
    plt.plot(tvec[1:], pwp_out['mld_exact'][1:], label='exact')
    plt.plot(tvec[1:], pwp_out['mld_exact2'][1:], label='exact2')
    plt.plot(tvec[1], pwp_out['mld'][1], 'ro', label='day0', ms=7)
    plt.plot(tvec[forcing['time']==1], pwp_out['mld'][forcing['time']==1], 'go', label='day1', ms=7)
    plt.plot(tvec[forcing['time']==5], pwp_out['mld'][forcing['time']==5], 'ko', label='day5', ms=7)
    plt.gca().invert_yaxis() 
    plt.xlabel('Time (days)')
    plt.ylabel('MLD (m)')
    plt.legend(loc=0)
    plt.grid(True)
    
    if 'dtime' in forcing:
        formatDates(plt.gca())
    
    if save_plots:     
        plt.savefig('plots/MLD_evolution_%s.pdf' %suffix, bbox_inches='tight')
        plt.close()
        
        
    #plot change in MLT, MLS and MLT elevation
    fig, axes = plt.subplots(3,1, sharex=True, figsize=(6.5, 8.5)) 
    axes[0].plot(tvec[1:], pwp_out['mlt'][1:]) 
    axes[0].set_ylabel('Temperature (C)')
    axes[0].set_title('Mixed Layer Temperature')
    axes[0].grid(True)
    
    axes[1].plot(tvec[1:], pwp_out['mls'][1:]) 
    axes[1].set_ylabel('Salinity (PSU)')
    axes[1].set_title('Mixed Layer Salinity')
    axes[1].grid(True)
    
    axes[2].semilogy(tvec[1:], pwp_out['mlt_elev'][1:]) 
    axes[2].set_ylabel('Temperature (C)')
    axes[2].set_xlabel('Time (days)')
    axes[2].set_title('MLT - T_fz')
    axes[2].grid(True) 
    
    if 'dtime' in forcing:
        formatDates(axes[2])
    
    if save_plots:     
        plt.savefig('plots/MLT_MLS_%s.pdf' %suffix, bbox_inches='tight')
        plt.close()
        
    
    #plot upper ocean temp, sal, oxy evolution 
    
    vble = ['temp', 'sal', 'ps']
    units = ['$^{\circ}$C', 'PSU', '$\mu$mol/kg']
    # clim = ([-1.5, 2.5], [34.0, 34.75], [200, 325])
    cmap = [plt.cm.coolwarm, plt.cm.RdYlGn_r, plt.cm.coolwarm]
    zlim = 250
    
    for i in range(3):
        fig = plt.figure(figsize=(7,9))
        ax1 = plt.subplot2grid((4,1), (0,0), colspan=1)
        ax2 = plt.subplot2grid((4,1), (1,0), rowspan=4)
    
        #plot Ice
        ax1.plot(tvec, 100*pwp_out['ice_thickness'], '-b')
        ax1.set_title('Ice thickness (cm)', fontsize=12)
        ax1.set_ylabel('Thickness (cm)', fontsize=12)
        ax1.get_xaxis().set_visible(False)
        ax1.grid(True)
    
        #plot ocean
        # im = ax2.pcolormesh(tvec, pwp_out['z'], pwp_out['temp'] -pwp_out['temp'][:,:1], vmin=-1.0, vmax=1.0, cmap=cmap[i])
        vble_anom = pwp_out[vble[i]]-pwp_out[vble[i]][:,:1]
        #ann_mean = pwp_out[vble[i]].mean(axis=1)
        #vble_anom = pwp_out[vble[i]]-ann_mean[:, np.newaxis]
        vmin, vmax =  -4*np.std(vble_anom), 4*np.std(vble_anom)
        #debug_here()
        im = ax2.pcolormesh(tvec, pwp_out['z'], vble_anom, vmin=vmin, vmax=vmax, cmap=cmap[i])
        ax2.plot(tvec, mld_exact2_sm, 'k')
        ax2.set_ylim(0,zlim)
        ax2.invert_yaxis() 
        ax2.set_ylabel('Depth (m)', fontsize=12)
        ax2.set_xlabel('Time (days)', fontsize=12)
        ax2.set_title('%s (%s) and MLD evolution' %(vble[i], units[i]), fontsize=12)
    
        cbar_ax = fig.add_axes([0.83, 0.10, 0.025, 0.58]) #make a new axes for the colorbar
        fig.subplots_adjust(right=0.8) #adjust sublot to make colorbar fit
        fig.subplots_adjust(hspace=0.3) 
        fig.colorbar(im, ax=ax2, cax=cbar_ax)
        
        if 'dtime' in forcing:
            formatDates(ax2)
    
        if save_plots:
            #pass
            plt.savefig('plots/MLD_%s_evolution_%s.png' %(vble[i], suffix), bbox_inches='tight')
            plt.close()
        
        
    # #plot salinity evolution
    # fig = plt.figure(figsize=(7,9))
    # ax1 = plt.subplot2grid((4,1), (0,0), colspan=1)
    # ax2 = plt.subplot2grid((4,1), (1,0), rowspan=4)
    #
    # #plot Ice
    # ax1.plot(tvec, 100*pwp_out['ice_thickness'], '-b')
    # ax1.set_title('Ice thickness (cm)', fontsize=12)
    # ax1.set_ylabel('Thickness (cm)', fontsize=12)
    # ax1.get_xaxis().set_visible(False)
    # ax1.grid(True)
    #
    # #plot ocean
    # im = ax2.pcolormesh(tvec, pwp_out['z'], pwp_out['sal']-pwp_out['sal'][:,:1], vmin=-0.5, vmax=0.5, cmap=plt.cm.RdYlGn_r)
    # ax2.plot(tvec, mld_exact2_sm, 'k')
    # # plt.plot(tvec, pwp_out['mld_approx'], 'w')
    # #ax2.set_ylim(0,250)
    # ax2.invert_yaxis()
    # ax2.set_ylabel('Depth (m)', fontsize=12)
    # ax2.set_xlabel('Time (days)', fontsize=12)
    # ax2.set_title('Salinity and MLD evolution', fontsize=12)
    #
    # cbar_ax = fig.add_axes([0.83, 0.10, 0.025, 0.58]) #make a new axes for the colorbar
    # fig.subplots_adjust(right=0.8) #adjust sublot to make colorbar fit
    # fig.subplots_adjust(hspace=0.3)
    # fig.colorbar(im, ax=ax2, cax=cbar_ax)
    #
    #
    # if save_plots:
    #     #pass
    #     plt.savefig('plots/MLD_sal_evolution_%s.png' %suffix, bbox_inches='tight')
    
    
    if showPlots:
        plt.show()
     

def custom_div_cmap(numcolors=11, name='custom_div_cmap', mincol='blue', midcol='white', maxcol='red'):
                    
    """ 
    Create a custom diverging colormap with three colors
    
    Default is blue to white to red with 11 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    """

    from matplotlib.colors import LinearSegmentedColormap 
    
    cmap = LinearSegmentedColormap.from_list(name=name, colors =[mincol, midcol, maxcol], N=numcolors)
    
    return cmap   
    
def save2nc(data_dict, fpath, dt_save=1, type='out'):
    
    #WARNING: This script is broken. Delete and re-do
    
    data_ds = xr.Dataset()
    data_dict_save = {}
    # zt_vars = ['temp', 'sal', 'uvel', 'vvel', 'dens']
    # t_vars = ['mld', 'F_atm', 'F_i', 'F_ocean_ice', 'ice_thickness', 'surf_ice_temp', 'mld_exact', 'mld_exact2']
    # z_vars = ['z']
    
    # if type=='out':
    #     zt_shape = data_dict['temp'].shape
    #     t_shape = data_dict['mld'].shape
    #     z_shape = data_dict['z'].shape
    # else:
    #     t_shape = data_dict['time'].shape
    
    tvec = data_dict['time']
    # len_tvec_save = int(np.floor(len(tvec)/dt_save))
    tvec_save = tvec[::dt_save]
    ti = range(0, len(tvec), dt_save)
    
    assert len(ti)==len(tvec_save), "Error. dimension mismatch."
    
    for key in data_dict:
        
        if isinstance(data_dict[key], (float, int)):
            data_ds[key] =  data_dict[key]
            data_dict_save[key] = data_dict[key]
            
        else:
            
            if type=='out':
                
                zt_shape = data_dict['temp'].shape
                t_shape = data_dict['mld'].shape
                z_shape = data_dict['z'].shape

                if data_dict[key].shape == t_shape:
                    dims = ('t',)
                    #need to take time mean of forcing. Can't subsample.
                    #data_mean = data_dict[key].reshape(-1, dt_save).mean(axis=1)
                    data_mean = []
                    for i in range(len(ti)):
                        #for when there is an incomplete day at the end
                        if ti[i]==ti[-1]:
                            data_mean.append(data_dict[key][ti[i]:].mean())             
                        else:
                            data_mean.append(data_dict[key][ti[i]:ti[i+1]].mean())
                            
                    data_mean = np.ma.array(data_mean)
                    data_ds[key] = (dims, data_mean)
                    data_dict_save[key] = np.array(data_mean)
    
                elif data_dict[key].shape == z_shape:
                    dims = ('z', )
                    data_ds[key] = (dims, data_dict[key])
                    data_dict_save[key] = data_dict[key]
    
                elif data_dict[key].shape == zt_shape:
                    dims = ('z', 't')                
                    # data_mean = np.zeros((len(data_dict['z']), len(tvec_save)))*np.nan
                    # for i in range(len(ti)):
                    #     #for when there is an incomplete day at the end
                    #     if ti[i]==ti[-1]:
                    #         data_mean[:,i] = data_dict[key][:, ti[i]:].mean(axis=1)
                    #     else:
                    #         data_mean[:,i] = data_dict[key][:, ti[i]:ti[i+1]].mean(axis=1)
                    #
                    # data_ds[key] = (dims, data_mean)
                    # data_dict_save[key] = data_mean

                    data_ds[key] = (dims, data_dict[key][:,::dt_save])
                    data_dict_save[key] = data_dict[key][:,::dt_save]

                else:
                    print("%s variable has unrecognized shape. Can't save to ncfile. Skipping..." %key)
                    continue
            else:
                
                if key=='absrb':
                    continue
                    
                dims = ('t',)
                
                data_mean = []
                for i in range(len(ti)):
                    #for when there is an incomplete day at the end
                    if ti[i]==ti[-1]:
                        data_mean.append(data_dict[key][ti[i]:].mean())             
                    else:
                        data_mean.append(data_dict[key][ti[i]:ti[i+1]].mean())
                        
                data_mean = np.ma.array(data_mean)
                data_ds[key] = (dims, data_mean)
                data_dict_save[key] = np.array(data_mean)

                # data_ds[key] = (dims, data_dict[key][::dt_save])
                # data_dict_save[key] = data_dict[key][::dt_save]

    data_dict['time'] = np.arange(len(tvec_save))
    data_ds['time'] = ('t', np.arange(len(tvec_save)))
    # debug_here()
    data_ds.to_netcdf(fpath)
    data_ds.close
    
    return data_dict_save
    


    
    
    
    
    
    
    