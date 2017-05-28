"""
File to contain ice model functions to go with PWP.py

reference: http://www.cesm.ucar.edu/models/atm-cam/docs/description/node37.html

"""

import numpy as np
import seawater as sw
import matplotlib.pyplot as plt
from IPython.core.debugger import Tracer
import timeit
import os
import PWP_helper as phf
import warnings

debug_here = Tracer()

#define constants
L_ice = 333e3 #Latent heat of fusion (J kg-1) - from Hyatt 2006
rho_ice = 920. #density of ice (kg/m3) - from Hyatt 2006
k_ice = 2. #thermal conductivity of sea ice (W m-1 K-1) - from Thorndike 1992
c_ice = 2e6/rho_ice #heat capacity of sea ice (J kg-1 K-1) - from Thorndike 1992/Hyatt 2006
c_sw = 4183.3 #heat capacity of seawater (J kg-1 K-1) - from Hyatt 2006
sal_ice = 4. #salinity of sea ice (PSU) - from Hyatt 2006
#alpha = 0.9 #ice fraction

melt_lyrs = 1 #TODO: play with this. PWP has issues with super thin, freshwater lenses.
bdry_lyr = 1
min_ice = 0 #m    #threshold for what is considered thin ice
T_fzi = 0.0 #freezing point of (pure) ice

#TODO: move these to params
sal_ref = 34.0
dens_ref = 1026.0
temp_fz = sw.fp(sal_ref, p=1)

#override_alpha = False


def iceGrowthModel_ode_v1(t, y, F_ai, F_oi, temp_fz, k_ice):
    
    """   
    # dT/dt =  (F_atm + F_ocean)/(c_ice*rho_ice*h)   (1)
    # dh/dt = -(k*(T-T_fz)/h - F_ocean)/(L*rho_ice)    (2)
    # where T is surface temp
        
    # in the model, let y1 = T and y2 = h   
    """
    
    dy0 = 2*(F_ai+F_oi)/(c_ice*y[1]*rho_ice)
    dy1 = (-k_ice*(y[0]-temp_fz)/y[1] - F_oi)/(L_ice*rho_ice)
    
    #debug_here()
    
    return [dy0, dy1]
    


def iceGrowthModel_ode_v2(t, y, temp_ice_surf, F_oi, temp_fz, k_ice):
    
    """   
    # dT/dt =  (F_atm + F_ocean)/(c_ice*rho_ice*h)   (1)
    # dh/dt = -(k*(T-T_fz)/h - F_ocean)/(L*rho_ice)    (2)
    # where T is surface temp
    
    In this model, the T is specified so the first eqn is eliminated.
    """
    
    F_i = -k_ice*(temp_ice_surf-temp_fz)/y[0]
    dh = (F_i - F_oi)/(L_ice*rho_ice)
    
    return [dh]

def iceGrowthModel_ode_v3(t, y, F_ai, F_oi): 
    
    """   
    Like version 2, but with surface temp chosen so that F_i == -F_ai 
    """
    F_i = -F_ai
    dh = (F_i - F_oi)/(L_ice*rho_ice)
    
    return [dh]
       

def get_ocean_ice_heat_flux(temp_sw, sal_sw, rho_sw, params):    
    
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    #temp_fz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    
    dT_surf = temp_sw[:bdry_lyr].mean() - temp_fz
    # F_sw = dT_surf*params['dz']*rho_sw0*c_sw*bdry_lyr #J per m2 per time step
    # F_sw_dt = F_sw/params['dt'] #in W/m2
    
    #mcphee turbulent heat flux parameterization
    c_d = 0.0056
    u_star = 0.01
    F_sw_dt = c_sw*dens_ref*c_d*u_star*dT_surf
    
    if F_sw_dt<0:
        #F_sw can turn negative if there is a large enough salt flux to lower layer1 temp past temp_fz(S=S_ref)
        #This is avoided if a constant reference salinity is used to compute temp_fz
        F_sw_dt=0.0
        
    

    return F_sw_dt

def create_initial_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, alpha, params):
    
    #temp_ice_surf_i not necessary I suppose
    if params['fix_alpha']:
        alpha = params['alpha']
        
    dz = params['dz']
    
    #temp_fz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    #rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    
    print("initiating ice growth...")
    ##create ice according to Hyatt 2006
    #first, create a thin layer of sea ice (eqn. 5.11, Hyatt 2006)
    h_ice_f = h_ice_i + dens_ref*c_sw*dz*(temp_fz-temp_sw[0])/(rho_ice*L_ice)
    
    #compute salinity change in top layer due to brine rejection (eqn. 5.12, Hyatt 2006)
    dh_ice = h_ice_f - h_ice_i
    dsal_sw = alpha*dh_ice*(sal_ref-sal_ice)/dz
    sal_sw[0] = sal_sw[0]+dsal_sw
    
    #set ice to freezing temp of water
    temp_ice_surf_f = temp_fz
    
    #set ocean surface to freezing temp
    temp_sw[0] = temp_fz
    
    
    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
    
    
def grow_ice_v1(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_ai, F_oi, alpha, params):
    
    """
    Function to evolve ice using when either F_ai>0 or when temp_ice_surf>T_fz
    
    """
    if params['fix_alpha']:
        #NOTE: alpha should be required in this script.
        alpha = params['alpha']
        
    F_aio = 0.0 #Leftover atm-ice heat flux that goes into warming the ocean.
        
    print("running basic ice algorithm...")
    
    #define surface density and freezing temp
    #temp_fz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    #rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
  
    #available heat: heat available to warm/melt of cool/grow ice:
    total_available_heat = (F_ai+F_oi)*params['dt']  #J/m2
    
    #warming_heat_sink: heat gain required to bring ice to freezing point (small compared to melting sink):
    if temp_ice_surf_i<temp_fz:
        warming_heat_sink = -(temp_ice_surf_i-temp_fz)*c_ice*h_ice_i #J/m2
    else:
        warming_heat_sink = 0.0
        
    #melting_heat_sink: heat gain required to completely melt ice, after it has warmed to the freezing point:
    melting_heat_sink = h_ice_i*rho_ice*L_ice #J/m2
    
    #total ice heat sink: warming_heat_sink + melting_heat_sink
    total_ice_heat_sink = warming_heat_sink + melting_heat_sink
    
    if total_available_heat>0:
    
        if total_available_heat<= warming_heat_sink:
            
            print("warming ice...")
            
            #use heat warm the ice. Assume linear profile.
            temp_ice_surf_f = temp_ice_surf_i + 2*total_available_heat/(c_ice*h_ice_i*rho_ice)
            #no change in ice thickness, since ice was not warmed to freezing temp
            dh_ice = 0
            
        else:
            
            """
            In this scenario total_available_heat exceeds heat needed to warm ice to freezing temp 
            so, use heat left over from warming to melt the ice and possibly warm the ocean.
            """
            
            print("melting ice...")
            
            available_heat_for_melt = total_available_heat-warming_heat_sink
            
            if available_heat_for_melt < melting_heat_sink:
                
                """
                In this scenario, there is not enough heat to completely melt the ice.
                """
                
                #compute change in ice thickness
                dh_ice = -available_heat_for_melt/(rho_ice*L_ice) #denom is J/m3. need to multiply by ice frac???
                
                #set ice surface temp to freezing temp
                temp_ice_surf_f = temp_fz
                
            else:
                
                """
                In this scenario, there is enough heat to warm the ice to its freezing temp, 
                AND completely melt it. That is, total_available_heat >= total_ice_heat_sink.
                
                
                NOTE:
                ice will completely melt in this iteration. Extra heat will warm ocean
                
                """
                
                print("Ice has completely melted.")
                dh_ice = -h_ice_i

                available_heat_for_ocean_warming = total_available_heat-total_ice_heat_sink
                assert available_heat_for_ocean_warming>=0, "Something went wrong. available_heat_for_ocean_warming should be sufficient to melt ice."
                
                #warm the ocean
                dT_surf = available_heat_for_ocean_warming/(params['dz']*dens_ref*c_sw)
                temp_sw[0] = temp_sw[0]+dT_surf
                
                #print the portion of the atm-ice heat flux that became an ocean warming flux.
                F_aio = available_heat_for_ocean_warming/params['dt']
                print("F_aio = %.2f" %F_aio)
                
                temp_ice_surf_f = np.nan #since there is no ice
                
                # TODO: delete soon
                # #first apply ocean heat flux
                # available_OCN_heat_for_ice_melt = F_oi*params['dt'] - warming_heat_sink
                # """
                # NOTE: available_OCN_heat_for_ice_melt will generally be positive. If it is negative,
                # it will add to the heat required to melt ice.
                # """
                #
                # heat_required_to_melt_ice = melting_heat_sink-available_OCN_heat_for_ice_melt
                # """
                # NOTE: This will generally be postive. If this is negative, then F_oi has completely
                # melted the ice. The leftover ocean heat for ice melt (i.e. excess ocean cooling) is
                # compensated by the available_ATM_heat_for_ocean_warming.
                # """
                #
                # #apply atmospheric heat flux
                # available_ATM_heat_for_ocean_warming = F_ai*params['dt'] - heat_required_to_melt_ice
                #
                # assert available_ATM_heat_for_ocean_warming>=0, "Something went wrong. F_ai should be sufficient to melt ice."
                #
                # #now warm the ocean
                # dT_surf = available_ATM_heat_for_ocean_warming/(params['dz']*dens_ref*c_sw)
                # temp_sw[0] = temp_sw[0]+dT_surf
                #
                # #the the portion of the atm-ice heat flux that became an ocean warming flux.
                # F_aio = available_ATM_heat_for_ocean_warming/params['dt']
                # print("F_aio = %.2f" %F_aio)
                #
                # temp_ice_surf_f = np.nan #since there is no ice    
                
                
        
    else:
        
        """
        if available heat is negative, grow ice according to Hyatt 2006 (eqn. 5.11). 
        """
        
        # message = "Oh oh. Unexpected error. Code should not come here..."
        # warnings.warn(message)
        # debug_here()

        print("growing ice following Hyatt 2006...")

        #compute ocean surface temp change due to heat loss through ice (combine with the above?)
        dT_surf = total_available_heat/(params['dz']*dens_ref*c_sw)
        temp_sw0 = temp_sw[0]+dT_surf

        #keep ice at the freezing point 
        temp_ice_surf_f = temp_fz

        #if there is enough heat to cool the ocean surface past its freezing point, use the remaining heat flux to create ice
        if temp_sw0 < temp_fz:
            dh_ice = dens_ref*c_sw*params['dz']*(temp_fz-temp_sw0)/(rho_ice*L_ice)
            temp_sw[0] = temp_fz
        else:
            #this is the case where the heat loss is not enough to cool the ocean past the freezing point.
            dh_ice = 0
            temp_sw[0] = temp_sw0

            
    #get final ice thickness
    h_ice_f = h_ice_i + dh_ice
    
    #debug_here()
  
    print("F_ai: %.2f, F_oi: %.2f, h_i=%.4f, dh: %.4f" %(F_ai, F_oi, h_ice_i, dh_ice))
    
    #debug_here()
    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw, F_aio
    
def grow_ice_v3(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_ai, F_oi, alpha, params):
    
    """
    This ice model uses specified surface ice temperatures.
    """
    from scipy.integrate import odeint
    import scipy.integrate as spi
    
    #derive/set additional parameters
    #rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    #temp_fz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    temp_ice_base = temp_fz
    F_i = 0.0 #heat flux through ice
    F_aio = 0.0 #heat flux from atmosphere left over after melting ice
    switch_algorithm = False

    if params['fix_alpha']:
        #use fixed alpha
        alpha = params['alpha']

    dz = params['dz']
    dt = params['dt']
    # alpha = params['alpha']
    
    assert h_ice_i >=0, "Error! negative ice thickness. Something went terribly wrong!"
    
    ### evolve sea ice ###
    if temp_ice_surf_i < temp_fz and F_ai<0: 
        print("grow ice with specified ice surface temperatures...")
    
        #define initial conditions and time step for ice growth model
        y_init = np.array([h_ice_i]) #initial ice thickness
        t_end = dt #end time (1 PWP time step)

        #load and intialize ice growth model
        dt_ice_model = 5 #seconds (sub-cycle PWP time step - necessary because heat flux through ice is inversely prop to h_i) 
        ode =  spi.ode(iceGrowthModel_ode_v2) #
        ode.set_f_params(temp_ice_surf_i, F_oi*dt_ice_model, temp_fz, k_ice*dt_ice_model)
    
        ode.set_integrator('lsoda')
        ode.set_initial_value(y_init, 0)
    
        ts = []
        ys = []
    
        #run ice growth model
        while ode.successful() and ode.t < t_end:
            ode.integrate(ode.t + dt_ice_model)
            h_ice = ode.y[0]
            ts.append(ode.t)
            ys.append(ode.y[0])
        
        
            #In the event F_oi large enough to completely melt ice, break out of this algorithm and use melt_ice algorithm
            ti = len(ys)-1
            if h_ice<=min_ice and ti>0 and ys[ti]<ys[ti-1]:
                print("ice is melting (has melted). Aborting grow_ice algorithm...")
                switch_algorithm=True

                #abort and switch to melt_ice algorithm.
                break


        if switch_algorithm:
            #this restarts the ice growth/melt process
            h_ice_f, temp_ice_surf_f, temp_sw, sal_sw, F_aio = grow_ice_v1(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_ai, F_oi, alpha, params)

        else:
            #get final ice thickness
            t = np.vstack(ts) #don't need this. Delete?
            y = np.vstack(ys)
            h_ice_f = y[-1,0]

            temp_ice_surf_f = temp_ice_surf_i
            
            
        t = np.vstack(ts) #don't need this. Delete?
        y = np.vstack(ys)
        h_ice_f = y[-1,0]

        temp_ice_surf_f = temp_ice_surf_i
            
        #get change in ice thickness
        dh_ice = h_ice_f - h_ice_i
    
        #find mean heat flux through ice
        F_i = -k_ice*(temp_ice_surf_f-temp_fz)/h_ice_f
        
        print("F_i: %.2f, F_oi: %.2f, h_i=%.2f, dh: %.4f" %(F_i, F_oi, h_ice_i, dh_ice))

    else:
        
        #use simple ice algorithm
        h_ice_f, temp_ice_surf_f, temp_sw, sal_sw, F_aio = grow_ice_v1(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_ai, F_oi, alpha, params)
        dh_ice = h_ice_f - h_ice_i
        
    
    assert h_ice_i >=0., "Error! negative ice thickness. Something went terribly wrong!"

    #cool ocean surface temp according to F_oi
    dT = F_oi*dt/(params['dz']*dens_ref*c_sw)
    temp_sw[0] = temp_sw[0]-dT
      
    
    #compute salinity change  
    dsal_sw = alpha*dh_ice*(sal_ref-sal_ice)/(dz*bdry_lyr)
    sal_sw[:bdry_lyr] = sal_sw[:bdry_lyr]+dsal_sw
    
    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw, F_i, F_aio

def ice_model_v4(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_ai, F_oi, alpha, params):
    
    """
    This ice model uses specificied surface fluxes and assumes the ice is in thermal equilibrium with the atmosphere
    
    UNDER CONSTRUCTION
    Do not use;
    """
    from scipy.integrate import odeint
    import scipy.integrate as spi
    
    #derive/set additional parameters
    #rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    #temp_fz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    temp_ice_base = temp_fz
    F_i = 0.0
    F_aio = 0.0
    switch_algorithm = False
    #temp_ice_surf_i = forcing['skt'] #experimenting...

    if params['fix_alpha']:
        alpha = params['alpha']

    dz = params['dz']
    dt = params['dt']
    # alpha = params['alpha']
    
    assert h_ice_i >=0., "Error! negative ice thickness. Something went terribly wrong!"
    
    
    #warming_sink: heat gain required to bring ice to freezing point:
    warming_heat_sink = -(temp_ice_surf_i-temp_fz)*c_ice*h_ice_i*alpha 
    #melting_heat_sink: heat gain required to completely melt ice of thickness, after it has warmed to the freezing point:
    melting_heat_sink = alpha*h_ice_i*rho_ice*L_ice
    #total ice heat sink: warming_heat_sink + melting_heat_sink
    total_ice_heat_sink = warming_heat_sink + melting_heat_sink
    
    
    ### Modify sea ice ###
    
    if abs(F_ai + F_oi) >= 0:
        
        #we split this up because, melting sea-ice is a bit tricky
        
        print("running ice model with specified ice-atmosphere fluxes..")
        h_ice_f = h_ice_i + dt*(-F_ai - F_oi)/(L_ice*rho_ice)  
        temp_ice_surf_f = F_ai*h_ice_f/k_ice+temp_fz 
        
        
        #get change in ice thickness
        dh_ice = h_ice_f - h_ice_i
    
        #debug_here()
        #find mean heat flux through ice
        F_i = -k_ice*(temp_ice_surf_f-temp_fz)/h_ice_f
        
        print("F_i: %.2f, F_oi: %.2f, h_i=%.2f, dh: %.4f" %(F_i, F_oi, h_ice_i, dh_ice))
        
    else:
        
        #use thin ice algorithm
        h_ice_f, temp_ice_surf_f, temp_sw, sal_sw, F_aio = grow_ice_v1(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_ai, F_oi, alpha, params)
        dh_ice = h_ice_f - h_ice_i

        
    
    assert h_ice_i >=0., "Error! negative ice thickness. Something went terribly wrong!"

    #cool ocean surface temp according to F_oi
    dT = F_oi*dt/(params['dz']*dens_ref*c_sw)
    temp_sw[0] = temp_sw[0]-dT
      
    
    #compute salinity change  
    dsal_sw = alpha*dh_ice*(sal_ref-sal_ice)/(dz*bdry_lyr)
    sal_sw[:bdry_lyr] = sal_sw[:bdry_lyr]+dsal_sw
    
    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw, F_i, F_aio        
    
        

        
    
    
    
    
    
    
    
    
    
    
    

    