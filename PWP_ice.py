"""
File to contain ice model functions to go with PWP.py

"""

import numpy as np
import seawater as sw
import matplotlib.pyplot as plt
from IPython.core.debugger import Tracer
import timeit
import os
import PWP_helper as phf

debug_here = Tracer()

#define constants
L_ice = 333e3 #Latent heat of fusion (J kg-1) - from Hyatt 2006
rho_ice = 920. #density of ice (kg/m3) - from Hyatt 2006
k_ice = 2. #thermal conductivity of sea ice (W m-1 K-1) - from Thorndike 1992
c_ice = 2e6/rho_ice #heat capacity of sea ice (J kg-1 K-1) - from Thorndike 1992/Hyatt 2006
c_sw = 4183. #heat capacity of seawater (J kg-1 K-1) - from Hyatt 2006
sal_ice = 5. #salinity of sea ice (PSU) - from Hyatt 2006

melt_lyrs = 5 #TODO: play with this. PWP has issues with super thin, freshwater lenses.
    

# def melt_all_ice(h_ice, temp_ice_surf, sal_sw, temp_sw):
#
#     print "Ice has completely melted..."
#
#     #compute freshening due to melt
#     dsal_sw_melt = h_ice*(sal_sw[:melt_lyrs].mean()-sal_ice)/(dz*melt_lyrs)
#     sal_sw[:melt_lyrs] = sal_sw[:melt_lyrs]+dsal_sw_melt
#
#     #compute temp loss due to melt
#     F_i = h_ice*rho_ice*L_ice #energy used to melt ice
#     dT_surf = F_i/(dz*rho_sw0*c_sw)
#     temp_sw[0] = temp_sw[0]-dT_surf
#
#     #reset ice
#     h_ice = 0.
#     temp_ice_surf = np.nan
#
#     #debug_here()
#     return h_ice, temp_ice_surf, temp_sw, sal_sw

def iceGrowthModel_ode(t, y, qnet, F_sw_dt, temp_swfz, k_ice):
    
    """   
    # dT/dt =  qnet/(c*rho*h)   (1)
    # dh/dt = -(k/hL)*(T-T_fz)   (2)
    # where T is surface temp
        
    # in the model, let y1 = T and y2 = h   
    """
    
    dy0 = 2*(qnet+F_sw_dt)/(c_ice*y[1]*rho_ice)
    dy1 = (-k_ice*(y[0]-temp_swfz)/y[1] - F_sw_dt)/L_ice
    
    return [dy0, dy1]
    

def get_ocean_ice_heat_flux(temp_sw, sal_sw, rho_sw, dz):    
    
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    
    dT_surf = temp_sw[0] - temp_swfz
    F_sw = dT_surf*dz*rho_sw0*c_sw #J per m2 per time step
    #F_sw_dt = F_sw/dt #in W/m2

    return F_sw
    
def melt_ice(h_ice_i, temp_ice_surf, temp_sw, sal_sw, rho_sw, q_net, dz, dt, source):
    
    """
    Function to melt ice
    
    source: 'ocean' or 'atm'
    
    q_net - heat flux (W/m2)
    
    only way ice temp changes here is if ALL the ice melts. Otherwise, the ice will be in 
    dis-equilibrium until next time step.
    
    """
    
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    
    if source == 'atm':
        q_net = q_net*dt
    
    
    #print "Applying ocean feedback..."  
    ice_melt_potential = q_net/(rho_ice*L_ice) #amount of ice that can be melted with available heat flux

    #check if ocean heat heat flux can completely melt current ice 
    if ice_melt_potential >= h_ice_i:
        print "%s heat flux has melted all ice..." %source
        #reset ice thickness and temperature
        temp_ice_surf = np.nan
        h_ice_f = 0.
        h_ice_melt = h_ice_i
        
    else:
        #ice survives only melt a portion of the ice with total ice melt potential
        h_ice_melt = ice_melt_potential
        h_ice_f = h_ice_i - h_ice_melt
    
    #compute freshening due to melt:
    dsal_sw_melt = h_ice_melt*(sal_sw[:melt_lyrs].mean()-sal_ice)/(dz*melt_lyrs)
    sal_sw[:melt_lyrs] = sal_sw[:melt_lyrs]+dsal_sw_melt
    
    #energy used to melt ice
    q_ice = h_ice_melt*rho_ice*L_ice 
    
    if source == 'ocean':
        #compute ocean surface temp loss due to melt:
        rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
        dT_surf = q_ice/(dz*rho_sw0*c_sw)
        temp_sw[0] = temp_sw[0]-dT_surf
    
        if ice_melt_potential < h_ice_i:
            assert np.abs(temp_sw[0]-temp_swfz)<0.001, "Something went wrong. Ocean surface not at freezing temp after incomplete ice melt."
            
        print "F_sw: %.2f, h_i: %.2f, dh: -%.4f" %(q_net/dt, h_ice_i, h_ice_melt)
        
    elif source == 'atm' and q_net>q_ice:
        
        #use excess atmospheric heat to warm ocean
        qnet_rem = q_net - q_ice
        dT_surf = qnet_rem/(dz*rho_sw0*c_sw)
        temp_sw[0] = temp_sw[0]-dT_surf
        
        print "F_atm: %.2f, h_i: %.2f, dh: %.4f" %(q_net/dt, h_ice_i, h_ice_melt)
    
    
    return h_ice_f, temp_ice_surf, temp_sw, sal_sw
        
   
def create_initial_ice(h_ice_i, temp_ice_surf, temp_sw, sal_sw, rho_sw, dz):
    
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    
    print "initiating ice growth..."
    ##create ice according to Hyatt 2006
    #first, create a thin layer of sea ice (eqn. 5.11, Hyatt 2006)
    h_ice_f = rho_sw0*c_sw*dz*(temp_swfz-temp_sw[0])/(rho_ice*L_ice)
    
    #compute salinity change in top layer due to brine rejection (eqn. 5.12, Hyatt 2006)
    dsal_sw = h_ice_f*(sal_sw[0]-sal_ice)/dz
    sal_sw[0] = sal_sw[0]+dsal_sw
    
    #set ice to freezing temp of water
    temp_ice_surf_f = temp_swfz
    
    #set ocean surface to freezing temp
    temp_sw[0] = temp_swfz
    
    
    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
    
    
            

def grow_existing_ice(temp_sw, sal_sw, rho_sw, z, dt, F_atm, temp_ice_surf_i, h_ice_i, params):
    """
    input:
    
    temp_sw - ocean temp profile
    sal_sw - ocean salinity profile
    rho_sw - ocean density profile
    z - vertical coordinates
    dt - time step in seconds
    F_atm - net surface heating (positive into the ocean)
    temp_ice_surf_i - initial surface ice temperature
    h_ice_i - initial ice thickness
    
    NOTE: We don't apply ocean heat flux here.
    """
    
    from scipy.integrate import odeint
    import scipy.integrate as spi
    
    #TODO: incorporate ice fraction
    
    #derive additional parameters
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    temp_ice_base = temp_swfz
    #F_atm_rem = 0. #heat left over after melting ice
    # h_ice_melt = 0. #ice melt due to ocean feedback
    F_sw = 0. #ocean heat flux
    dz = z[1]-z[0]
    
    assert h_ice_i >=0., "Error! negative ice thickness. Something went terribly wrong!"

    if F_atm > 0: 
        #if qnet is positive melt ice, the straightforward way.
        h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = melt_ice(h_ice_i, temp_ice_surf, temp_sw, sal_sw, rho_sw, F_atm, dz, dt, source='atm')
        
    else:
        #define initial conditions and time step for ice growtg model
        y_init = np.array([temp_ice_surf_i, h_ice_i])
        t_end = dt  
        t_ice_model = np.linspace(0, t_end, 1e5)
        dt_ice_model = t_ice_model[1]-t_ice_model[0]
    
        #load and intialize ice growth model
        ode =  spi.ode(iceGrowthModel_ode)
        #ode.set_integrator('vode', nsteps=500, method='adams') # BDF method suited to stiff systems of ODEs
        ode.set_integrator('lsoda')
        ode.set_initial_value(y_init, 0)
        ode.set_f_params(F_atm*dt_ice_model, F_sw*dt_ice_model, temp_swfz, k_ice*dt_ice_model)
        ts = []
        ys = []
        
        #debug_here()

        #run ice growth model
        while ode.successful() and ode.t < t_end:
            ode.integrate(ode.t + dt_ice_model)
            ts.append(ode.t)
            ys.append(ode.y)

        # while ode.successful() and ode.t < t_end: ode.integrate(ode.t + tstep_ice); ts.append(ode.t); ys.append(ode.y)
        t = np.vstack(ts)
        y = np.vstack(ys)
        temp_ice_surf_f = y[-1,0] 
        h_ice_f = y[-1,1]
        
        #compute ice temp and thickness change
        dh_ice = h_ice_f - h_ice_i
        dtemp_ice_surf = temp_ice_surf_f - temp_ice_surf_i
    
        print "F_atm: %.2f, h_i=%.2f, dh: %.4f" %(F_atm, h_ice_i, dh_ice)
    
        #check if all ice has melted
        assert h_ice_f >= h_ice_i, "Something weird happened. Ice seemed to have shrunk with negative fluxes???"
        # if h_ice_f <= h_ice_i:
        #
        #     print "Warning: ice model produced negative ice thickness. Not good..."
        #     print "Reseting to zero ice..."
        #
        #     h_ice, temp_ice_surf, temp_sw, sal_sw = melt_all_ice(h_ice, temp_ice_surf, sal_sw, temp_sw)
        #
        #     return h_ice, temp_ice_surf, temp_sw, sal_sw, qnet_rem
        
    
        #compute salinity change in top layer due to brine rejection (or ice melt?) (eqn. 5.12, Hyatt 2006)
        num_lyrs = 5 #TODO: play with this. PWP has issues with super thin, freshwater lenses.
        dsal_sw = dh_ice*(sal_sw[:num_lyrs].mean()-sal_ice)/(dz*num_lyrs)
        sal_sw[:num_lyrs] = sal_sw[:num_lyrs]+dsal_sw

        h_ice = h_ice_f
        temp_ice_surf = temp_ice_surf_f
        
    return h_ice, temp_ice_surf, temp_sw, sal_sw
    
    
    # if h_ice_i == 0.:
    #     print "initiating ice growth..."
    #     ##create ice according to Hyatt 2006
    #     #first, create a thin layer of sea ice (eqn. 5.11, Hyatt 2006)
    #     h_ice_f = rho_sw0*c_sw*dz*(temp_swfz-temp_sw[0])/(rho_ice*L_ice)
    #
    #     #compute salinity change in top layer due to brine rejection (eqn. 5.12, Hyatt 2006)
    #     dsal_sw = h_ice_f*(sal_sw[0]-sal_ice)/dz
    #     sal_sw[0] = sal_sw[0]+dsal_sw
    #
    #     #set ice to freezing temp of water
    #     temp_ice_surf_f = temp_swfz
    #
    #     #set ocean surface to freezing temp
    #     temp_sw[0] = temp_swfz

    #elif h_ice_i > 0:
        
        
           
        # h_ice_i = h_ice #initial ice thickness
        # temp_ice_surf_i = temp_ice_surf #initial surface ice temperature
        
        
        # if ocean_fb_ON:
        #     # Find ocean sensible heat flux.
        #     # That is, the heat flux required to bring surface temp to freezing.
        #     dT_surf = temp_sw[0] - temp_swfz
        #     F_sw = dT_surf*dz*rho_sw0*c_sw
        #     F_sw_dt = F_sw/dt
        #

        #check if ocean heat heat flux can completely melt current ice
        # h_ice_melt = F_sw/(rho_ice*L_ice) #amount of ice that can be melted
        # if h_ice_melt >= h_ice:
        #
        #     h_ice, temp_ice_surf, temp_sw, sal_sw =  melt_all_ice(h_ice, temp_ice_surf, sal_sw, temp_sw)
        #     return h_ice, temp_ice_surf, temp_sw, sal_sw, qnet_rem
            
        # else:
        #
        #     temp_sw[0] = temp_swfz

                # print "Ocean mixing has melted ice..."
                #
                # #compute freshening due to melt
                # # dsal_sw_melt = h_ice*(sal_sw[0]-sal_ice)/dz
                # # sal_sw[0] = sal_sw[0]-dsal_sw_melt
                # #
                #
                # dsal_sw_melt = h_ice*(sal_sw[:melt_lyrs].mean()-sal_ice)/(dz*melt_lyrs)
                # sal_sw[:melt_lyrs] = sal_sw[:melt_lyrs]+dsal_sw_melt
                #
                # #compute temp loss due to melt
                # F_i = h_ice*rho_ice*L_ice #energy used to melt ice
                # dT_surf = F_i/(dz*rho_sw0*c_sw)
                # temp_sw[0] = temp_sw[0]-dT_surf
                #
                # #reset ice
                # h_ice = 0.
                # temp_ice_surf = np.nan

                #debug_here()
            #     return h_ice, temp_ice_surf, temp_sw, sal_sw, qnet_rem
            #
            # else:
            #
            #     temp_sw[0] = temp_swfz
                

        #debug_here()
        #h_ice_i = h_ice #initial ice thickness after ocean melt
        

        


    

        

    
        
        
        
        # plt.figure()
        # plt.subplot(211)
        # plt.plot(y[:,0]); plt.ylabel('Ice temperature (C)')
        # plt.subplot(212)
        # plt.plot(y[:,1]); plt.ylabel('Ice thickness (m)')
        
         
    
    

    