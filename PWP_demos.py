import PWP
import PWP_helper as phf


def demo1():
    """
    Example script of how to run the PWP model. 
    This run uses summertime data from the Beaufort gyre
    """
    
    forcing_fname = 'beaufort_met.nc'
    prof_fname = 'beaufort_profile.nc' 
    print("Running Test Case 1 with data from Beaufort gyre...")
    
    
    #define param modifications
    p={}
    p['dens_option'] = 'pdens'
    p['plot_zlim'] = 100 #meters
    
    suffix='demo1_nodiff' #string to append to all plots names

    forcing, pwp_out = phf.run_PWP(met_data=forcing_fname, prof_data=prof_fname, makeLivePlots=False, suffix=suffix, save_plots=True, param_mods=p)
    
    return forcing, pwp_out, suffix

def demo2(dopt='pdens', winds_ON=True, emp_ON=True):
    
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
    p['dt_hr'] = 6
    p['dens_option'] = dopt
    p['examine_stabilized_plot'] = False
    p['quiet_mode'] = False
    p['plots2make'] = [1, 6]
    
    #surface forcing controls
    p['winds_ON'] = winds_ON 
    p['emp_ON'] = emp_ON
    
    #create file name based on param settings
    if p['winds_ON']:
        wind_str = ''
    else:
        wind_str = '_noWINDS'
        
    if p['emp_ON']:
        emp_str = ''
    else:
        emp_str = '_noEMP'
    suffix = 'demo2_1e6diff%s%s'%(wind_str, emp_str)

    forcing, pwp_out = phf.run_PWP(met_data=forcing_fname, prof_data=prof_fname, suffix=suffix, save_plots=True, param_mods=p)
    
    return forcing, pwp_out, suffix 

def demo3(iceMod=1):
    
    """
    Example run that involves sea-ice
    This run is initialized with an early winter profile from Maud Rise and uses NCEP fluxes. 
    """
    
    p={}
    p['rkz']=1e-5
    p['dz'] = 1
    p['dt_hr'] = 6
    p['max_depth'] = 1000
    p['ice_ON'] = True
    p['winds_ON'] = True
    p['emp_ON'] = False
    p['alpha'] = 0.95
    p['dens_option'] = 'pdens'
    # p['fix_alpha'] = True
    p['mld_thresh'] = 0.0001
    p['use_Bulk_Formula'] = True
    p['qnet_offset'] = 0 #W/m2 (default is zero)
    p['iceMod'] = iceMod #0: use ice_model_0(), 1: ice_model_T()
    p['gradMix_ON'] = False
    p['examine_stabilized_plot'] = False #generally you want this to be ON
    p['quiet_mode'] = False #to suppress out
    p['plots2make'] = [2, 6, 10, 11, 12] #plots to make 
    
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
        iceMod_str = '_iceModT'
    elif p['iceMod']==0:
        iceMod_str = '_iceMod0'
    elif p['iceMod']==2:
        iceMod_str = '_iceModF'
    
        
    
    if p['gradMix_ON']:
        gradMix_str = ''
    else:
        gradMix_str = '_noGradMix'
        
    
    fnum = 9094#'0068' #9099
    p1 = 10#25 #10
    p2 = 25#28 #12
    nump = p2-p1
    met_data = 'NCEP_forcing_for_f%s_p%s-%s.nc' %(fnum, p1, p2)
    prof_data = 'float%s_%s_%s.nc' %(fnum, p1, p2)
    suffix = '%s_%sc%s%s%s%s_alpha%s%s%s%s' %(fnum, nump, ice_str, wind_str, emp_str, qflux_str, p['alpha'], q_offset_str, iceMod_str, gradMix_str)
    
    
    forcing, pwp_out = phf.run_PWP(met_data=met_data, prof_data=prof_data, param_mods=p, suffix=suffix, save_plots=True)
    
    return forcing, pwp_out, suffix
    
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
    p['dt'] = 6
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
    
    forcing, pwp_out = phf.run_PWP(met_data=met_data, prof_data=prof_data, param_mods=p, suffix=suffix, save_plots=True)
    
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
    p['dt'] = 6
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
    
    forcing, pwp_out = phf.run_PWP(met_data=met_data, prof_data=prof_data, param_mods=p, suffix=suffix, save_plots=True)
    
    return forcing, pwp_out, suffix    
  
    