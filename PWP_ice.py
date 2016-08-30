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
sal_ice = 4. #salinity of sea ice (PSU) - from Hyatt 2006
#alpha = 0.9 #ice fraction

melt_lyrs = 2 #TODO: play with this. PWP has issues with super thin, freshwater lenses.
bdry_lyr = 1
thin_ice = 0.001 #m    
T_fzi = 0.0 #freezing point of (pure) ice
override_alpha = True


# def melt_all_ice(h_ice, temp_ice_surf, sal_sw, temp_sw):
#p[]
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

def iceGrowthModel_ode(t, y, F_atm, F_ocean, temp_swfz, k_ice):
    
    """   
    # dT/dt =  (F_atm + F_ocean)/(c_ice*rho_ice*h)   (1)
    # dh/dt = -(k*(T-T_fz)/h - F_ocean)/(L*rho_ice)    (2)
    # where T is surface temp
        
    # in the model, let y1 = T and y2 = h   
    """
    
    dy0 = 2*(F_atm+F_ocean)/(c_ice*y[1]*rho_ice)
    dy1 = (-k_ice*(y[0]-temp_swfz)/y[1] - F_ocean)/(L_ice*rho_ice)
    
    return [dy0, dy1]
    

def iceGrowthModel_ode_v2(t, y, temp_ice_surf, F_ocean, temp_swfz, k_ice):
    
    """   
    In this version, the surface temperature is specified.
    """
    
    F_i = -k_ice*(temp_ice_surf-temp_swfz)/y[0]
    dy0 = (F_i - F_ocean)/(L_ice*rho_ice)
    
    return [dy0]
    

def get_ocean_ice_heat_flux(temp_sw, sal_sw, rho_sw, params):    
    
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    
    dT_surf = temp_sw[:bdry_lyr].mean() - temp_swfz
    F_sw = dT_surf*params['dz']*rho_sw0*c_sw*bdry_lyr #J per m2 per time step
    if F_sw<0:
        #F_sw can turn negative if there is freshening event that increases the freezing point of water
        F_sw=0.0
        
    F_sw_dt = F_sw/params['dt'] #in W/m2

    return F_sw_dt
    
def melt_ice(h_ice_i, temp_ice_surf, temp_sw, sal_sw, rho_sw, q_net, dz, dt, source):
    #TODO: delete this, if no longer needed
    """
    Defunct.
    
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

    #check if heat flux can completely melt current ice 
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
        temp_sw[0] = temp_sw[0]-dT_surf #shouldn't this be positive temp_sw[0]+dT_surf ???
        
        print "F_atm: %.2f, h_i: %.2f, dh: -%.4f" %(q_net/dt, h_ice_i, h_ice_melt)
    
    
    return h_ice_f, temp_ice_surf, temp_sw, sal_sw
        
   
def create_initial_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, params):
    
    #temp_ice_surf_i not necessary I suppose
    if override_alpha:
        alpha = params['alpha']
        
    dz = params['dz']
    
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
    
    
def modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, alpha, params):
    
    """
    Use this method to:
    
    1. modify thin ice up to some threshold. 
    2. melt thick or thin ice if available heat flux is sufficient. ODE solver cannot handle this.
    
    """
    if override_alpha:
        alpha = params['alpha']
    
    #define surface density and freezing temp
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)

    
    #available heat: heat available to warm/melt of cool/grow ice:
    total_available_heat = alpha*(F_atm+F_ocean)*params['dt']  
    #warming_sink: heat gain required to bring ice to freezing point:
    warming_heat_sink = (temp_ice_surf_i-temp_swfz)*c_ice*h_ice_i*alpha 
    #melting_heat_sink: heat gain required to completely melt ice of thickness, after it has warmed to the freezing point:
    melting_heat_sink = alpha*h_ice_i*rho_ice*L_ice
    #total ice heat sink: warming_heat_sink + melting_heat_sink
    total_ice_heat_sink = warming_heat_sink + melting_heat_sink
    
    if total_available_heat>0:
        
        #apply net ocean heat flux through leads
        # dT_surf1 = (1-alpha)*(F_atm)*params['dt']/(params['dz']*rho_sw0*c_sw) #TODO: should ocean heat fluxes apply here??
        # temp_sw[0] = temp_sw[0]+dT_surf1
        
        
        if total_available_heat<= warming_heat_sink:
            
            print "warming ice..."
            
            #use heat warm the ice
            temp_ice_surf_f = temp_ice_surf_i + 2*total_available_heat/(alpha*c_ice*h_ice_i*rho_ice)
            #no change in ice thickness, since ice was not warmed to freezing temp
            dh_ice = 0
            #reset ocean surface to freezing temp:
            temp_sw[:bdry_lyr] = temp_swfz
            
        else:
            
            """
            In this scenario total_available_heat exceeds heat needs to warm ice to freezing temp 
            Hence, use heat left over from warming to melt the ice and possibly warm the ocean.
            """
            
            print "melting ice..."
            
            available_heat_for_melt = total_available_heat-warming_heat_sink
            
            if available_heat_for_melt <= melting_heat_sink:
                
                """
                In this scenario, there is enough heat to warm the ice to its freezing temp, 
                but not enough to completely melt it.
                """
                
                #compute change in ice thickness
                dh_ice = -available_heat_for_melt/(rho_ice*L_ice*alpha) #denom is J/m3. need to multiply by ice frac
                #reset ocean surface temp:
                temp_sw[:bdry_lyr] = temp_swfz
                #set ice surface temp to freezing
                temp_ice_surf_f = temp_swfz
                
            else:
                
                """
                In this scenario, there is enough heat to warm the ice to its freezing temp, 
                AND completely melt it. That is, total_available_heat > total_ice_heat_sink
                """
                
                #first find available ATMOSPHERIC heat.
                #here we multiply by alpha because we already applied (1-alpha)*(F_atm+F_ocean) to the ocean surface
                available_ATM_heat_for_ocean_warming = alpha*F_atm*params['dt'] - (total_ice_heat_sink-alpha*F_ocean*params['dt'])
                
                #completely melt ice and use leftover ATMOSPHERIC heat to warm the ocean
                dh_ice = -h_ice_i
                dT_surf = available_ATM_heat_for_ocean_warming/(params['dz']*rho_sw0*c_sw)
                temp_sw[0] = temp_sw[0]+dT_surf
                
                #set ice surface temp to freezing
                temp_ice_surf_f = np.nan #since there is no ice
                
                print "Ice has completely melted."
                
        
    elif total_available_heat <=0 and alpha>0:
        
        """
        if available heat is negative, grow ice according to Hyatt 2006 (eqn. 5.11). 
        
        With this approach, we assume that ice is essentially transparent to radiative
        heat fluxes. The intent is to mitigate the rapid ice growth that occurs when 
        we apply large negative heat fluxes to very thin ice.
        
        With this model, ice continues to grow/accumulate until it exceeds the thin ice
        threshold. 
        
        """
        
        print "growing thin ice..."
        
        #set ocean to freezing point. If sea ice exists, ocean surface needs to be at freezing temp
        temp_sw[0] = temp_swfz 
        #if total_available_heat < 0, that means there is enough atm cooling to bring ocean temp to the freezing point
        #what's leftover, we use to creat sea-ice
        
        #apply net heat flux through leads
        dT_surf1 = (1-alpha)*(F_atm)*params['dt']/(params['dz']*rho_sw0*c_sw)
        temp_sw0 = temp_sw[0]+dT_surf1
    
        #compute ocean surface temp change due to heat loss through ice (combine with the above?)
        dT_surf = total_available_heat/(params['dz']*rho_sw0*c_sw)
        temp_sw0 = temp_sw0+dT_surf
        
        #keep ice at the freezing point (might be redundant)
        temp_ice_surf_f = temp_swfz
        
        #if there is enough heat to cool the ocean surface past its freezing point, use the remaining heat flux to create ice
        if temp_sw0 < temp_swfz:       
            dh_ice = rho_sw0*c_sw*params['dz']*(temp_swfz-temp_sw0)/(rho_ice*L_ice)
            
        else:
            #this is the case where the heat loss is not enough to cool the ocean past the freezing point.
            dh_ice = 0
            
            
    else:
        
        assert alpha==0.0, "non-zero. ice fraction."
        
        temp_ice_surf_f = np.nan
        
        #apply net heat flux to open ocean
        dT_surf = F_atm/(params['dz']*rho_sw0*c_sw)
        temp_sw0 = temp_sw[0]+dT_surf
        
        #if there is enough heat to cool the ocean surface past its freezing point, use the remaining heat flux to create ice
        if temp_sw0 < temp_swfz:       
            dh_ice = rho_sw0*c_sw*params['dz']*(temp_swfz-temp_sw0)/(rho_ice*L_ice)
            temp_sw[0] = temp_swfz
        else:
            #this is the case where the heat flux is not enough to cool the ocean past the freezing point.
            temp_sw[0] = temp_sw0
            dh_ice = 0
        
        
            
    #get final ice thickness
    h_ice_f = h_ice_i + dh_ice
    
    #debug_here()
  
    #compute salinity change in top layer due to brine rejection or ice melt (eqn. 5.12, Hyatt 2006)
    #NOTE: alpha>0 should not affect dh_ice directly. They should cancel out in ice equations above. 
    #Need it here to ensure correct surface salinity change.
    if dh_ice<0:
        dsal_sw = alpha*dh_ice*(sal_sw[0]-sal_ice)/params['dz'] 
    else:
        #unlike the case above, here we assume that when thin ice grow, it completely fills the grid box
        #NOTE: this assumption is a bit sketchy.
        dsal_sw = dh_ice*(sal_sw[0]-sal_ice)/params['dz'] 
        
    sal_sw[0] = sal_sw[0]+dsal_sw
    
    print "F_atm: %.2f, F_ocean: %.2f, h_i=%.4f, dh: %.4f" %(F_atm, F_ocean, h_ice_i, dh_ice)
    
    #debug_here()
    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
    
     

def modify_existing_ice(temp_ice_surf_i, h_ice_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, alpha, skt, params):
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
    alpha - ice fraction
    """
    
    from scipy.integrate import odeint
    import scipy.integrate as spi
    
    #TODO: current implementation of alpha is problematic. 
    #Could have situations where alpha=0 but model belives there should be ice. Could this be a side effect of advection?
    
    limit_ice_temp = True
    print_ice_temp_warning = True
    switch_algorithm = False
    
    #derive/set additional parameters
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    temp_ice_base = temp_swfz
    #F_atm_rem = 0. #heat left over after melting ice
    # h_ice_melt = 0. #ice melt due to ocean feedback
    # F_sw = 0. #ocean heat flux
    dz = params['dz']
    dt = params['dt']
    
    if override_alpha:
        alpha = params['alpha']
    
    assert h_ice_i >=0., "Error! negative ice thickness. Something went terribly wrong!"
    
    if alpha==0:
        
        print "ice-frac = 0. reseting ice thickness and ice surface temp..."
        h_ice_i = 0.0
        temp_ice_surf_i = np.nan
        # h_ice_f, temp_ice_surf_f, temp_sw, sal_sw  = create_initial_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, params)
        #
        # return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw


    if h_ice_i <= thin_ice:
        #modeling the conductive heat flux through thin ice is numerically unstable
        #use simplified approach that treats ice as transparent, except when positive heating is applied.
        h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, alpha, params)
        
        return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
        
        
    else:
        
        #available heat: heat required to warm/melt or cool/grow ice:
        available_heat = alpha*(F_atm+F_ocean)*params['dt']  
        #warming_sink: heat gain required to bring ice to freezing point:
        warming_heat_sink = (temp_ice_surf_i-temp_swfz)*c_ice*h_ice_i*alpha
        #melting_heat_sink: heat gain required to completely* melt ice of thickness, after it has warmed to the freezing point:
        melting_heat_sink = (h_ice_i-thin_ice)*rho_ice*L_ice*alpha #*or melt ice to the point where it becomes thin
        #total ice heat sink: warming_heat_sink + melting_heat_sink
        total_ice_heat_sink = warming_heat_sink + melting_heat_sink
        
        
        
        if available_heat >= total_ice_heat_sink:
        
            # if there is enough avaiable heat to create thin ice or completely melt the ice,
            # use thin ice algorithm. This is to circumvent the numerical issues that arise with very thin ice.
            # (this if-block is separated from the previous for the sake of clarity.)
        
            h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, alpha, params)
        
            return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
            
        else:
        
            print "running thick ice algorithm..."
        
            #if the ice is "thick", we can treat is using the Thorndike algorithm
            
            # if available_heat < total_ice_heat_sink, then ocean temp needs to be set to the freezing point
            # moreover, temp_sw[0] needs to be set to zero if ice is present
            temp_sw[0] = temp_swfz
        
            #first apply heat flux through leads
            dT_surf1 = (1-alpha)*(F_atm+F_ocean)/(params['dz']*rho_sw0*c_sw)
            
            #temp_sw[0] = temp_sw[0]+dT_surf1 #TODO: what if this causes the ocean temp to go below freezing? FIXED
            #first apply heat flux through leads
            if dT_surf1 < 0:
                #grow ice
                dh_ice_leads = rho_sw0*c_sw*params['dz']*(dT_surf1)/(rho_ice*L_ice)
                h_ice_i = h_ice_i + dh_ice_leads
                #create brine
                dsal_sw_leads = alpha*dh_ice_leads*(sal_sw[0]-sal_ice)/params['dz'] 
                sal_sw[0] = sal_sw[0]+dsal_sw_leads
                #debug_here()
            else:
                #warm ocean
                temp_sw[0] = temp_sw[0]+dT_surf1
                
            
            #define initial conditions and time step for ice growth model
            y_init = np.array([temp_ice_surf_i, h_ice_i])
            t_end = dt  
            t_ice_model = np.linspace(0, t_end, 1e4)
            dt_ice_model = t_ice_model[1]-t_ice_model[0]

            #load and intialize ice growth model
            ode =  spi.ode(iceGrowthModel_ode)
            #ode.set_integrator('vode', nsteps=500, method='adams') # BDF method suited to stiff systems of ODEs
            ode.set_integrator('lsoda')
            ode.set_initial_value(y_init, 0)
            ode.set_f_params(F_atm*dt_ice_model*alpha, F_ocean*dt_ice_model*alpha, temp_swfz, k_ice*dt_ice_model)
            ts = []
            ys = []
    
            #debug_here()

            #run ice growth model
            while ode.successful() and ode.t < t_end:
                ode.integrate(ode.t + dt_ice_model)
                temp_ice = ode.y[0] 
                h_ice = ode.y[1]
            
                if temp_ice >=temp_swfz or h_ice<=thin_ice:
                    print "ice has become too thin or warm. May become numerically unstable. Switching to thin ice algorithm..."
                    switch_algorithm=True

                    #abort and switch to thin ice algorithm.??
                    break
                    
                else:
                    
                    if limit_ice_temp:
                        if temp_ice < skt:
                            ode.y[0] = skt
                            
                            if print_ice_temp_warning:
                                print "enforcing lower limit on ice temp"
                                print_ice_temp_warning = False
                    
                    ts.append(ode.t)
                    ys.append(ode.y)


        

            if switch_algorithm:
            
                #use simplified approach that treats ice as transparent, except when positive heating is applied.
                h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, alpha, params)
        
                return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
            
            
            else:
                t = np.vstack(ts)
                y = np.vstack(ys)
                temp_ice_surf_f = y[-1,0] 
                h_ice_f = y[-1,1]
            
                #set surface temp to freezing
                temp_sw[0] = temp_swfz

            #compute ice temp and thickness change
            dh_ice = h_ice_f - h_ice_i
            dtemp_ice_surf = temp_ice_surf_f - temp_ice_surf_i

            print "F_atm: %.2f, F_ocean: %.2f, h_i=%.2f, dh: %.4f" %(F_atm, F_ocean, h_ice_i, dh_ice)

            #check if all ice has melted
            assert h_ice_f >= 0, "Something weird has happened. Negative ice thickness..."
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
            dsal_sw = alpha*dh_ice*(sal_sw[:num_lyrs].mean()-sal_ice)/(dz*num_lyrs)
            sal_sw[:num_lyrs] = sal_sw[:num_lyrs]+dsal_sw

        
            return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw
    
def ice_model_v3(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, params, forcing):
    
    """
    This ice model uses specified surface ice temperatures.
    
    UNDER CONSTRUCTION
    """
    from scipy.integrate import odeint
    import scipy.integrate as spi
    
    #derive/set additional parameters
    rho_sw0 = rho_sw[0] #density of seawater in model layer 1 (kg/m3)
    temp_swfz = sw.fp(sal_sw[0], p=1) #freezing point of seawater at given salinity
    temp_ice_base = temp_swfz
    F_i = np.nan
    switch_algorithm = False
    #temp_ice_surf_i = forcing['skt'] #experimenting...


    dz = params['dz']
    dt = params['dt']
    alpha = params['alpha']
    
    assert h_ice_i >=0., "Error! negative ice thickness. Something went terribly wrong!"
    
    
    #warming_sink: heat gain required to bring ice to freezing point:
    warming_heat_sink = (temp_ice_surf_i-temp_swfz)*c_ice*h_ice_i*alpha 
    #melting_heat_sink: heat gain required to completely melt ice of thickness, after it has warmed to the freezing point:
    melting_heat_sink = alpha*h_ice_i*rho_ice*L_ice
    #total ice heat sink: warming_heat_sink + melting_heat_sink
    total_ice_heat_sink = warming_heat_sink + melting_heat_sink
    
    
    ### first apply heat flux through leads ###
    # dT_surf = (1-alpha)*(F_atm*params['dt'])/(params['dz']*rho_sw0*c_sw) #TODO: Should this include F_ocean? No.
    # temp_sw0 = temp_sw[0]+dT_surf
    # #TODO: apply this to the entire ML???
    # #TODO: need to fix
    # if temp_sw0 < temp_swfz:
    #
    #     #set surf temp to freezing
    #     temp_sw[0] = temp_swfz
    #
    #     #grow ice
    #     dh_ice_leads = rho_sw0*c_sw*params['dz']*(temp_swfz-temp_sw0)/(rho_ice*L_ice)
    #     h_ice_i = h_ice_i + dh_ice_leads
    #
    #     #create brine
    #     dsal_sw_leads = alpha*dh_ice_leads*(sal_sw[0]-sal_ice)/params['dz']
    #     sal_sw[0] = sal_sw[0]+dsal_sw_leads
    #
    #     debug_here()
    #     #need to adjust F_ocean
    #     #F_ocean = F_ocean - (1-alpha)*(F_atm) #?????
    #
    # else:
    #     #cool/warm ocean.
    #     temp_sw[0] = temp_sw[0]+dT_surf
        
    #debug_here()  
    ### Modify sea ice ###
    #TODO: eventually this should work for all temperatures
    if temp_ice_surf_i < -2 and F_atm<0:
        
        #print "running sea ice algorithm..."
        print "running ice model with specified surface temperatures..."
        
        #define initial conditions and time step for ice growth model
        y_init = np.array([h_ice_i])
        t_end = dt  
        t_ice_model = np.linspace(0, t_end, 1e4)
        dt_ice_model = t_ice_model[1]-t_ice_model[0]

        #load and intialize ice growth model
        ode =  spi.ode(iceGrowthModel_ode_v2)
        ode.set_integrator('lsoda')
        ode.set_initial_value(y_init, 0)
        ode.set_f_params(temp_ice_surf_i, F_ocean*dt_ice_model*alpha, temp_swfz, k_ice*dt_ice_model)
        ts = []
        ys = []
        
        #run ice growth model
        while ode.successful() and ode.t < t_end:
            ode.integrate(ode.t + dt_ice_model)
            h_ice = ode.y[0]
            ts.append(ode.t)
            ys.append(ode.y[0])
            
            
            #TODO: Add a clause that causes the loop to break if ice is thin AND melting
            #debug_here()
            ti = len(ys)-1
            if h_ice<=thin_ice and ti>0 and ys[ti]<ys[ti-1]:
                print "ice is melting and has become too thin. Switching to thin ice algorithm..."
                switch_algorithm=True

                #abort and switch to thin ice algorithm.??
                break

        
        if switch_algorithm:
            #this restarts the ice growth/melt process
            h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, alpha, params)
            
            return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw, F_i
            
        else:
            t = np.vstack(ts)
            y = np.vstack(ys)
            h_ice_f = y[-1,0]
        
        #get change in ice thickness
        dh_ice = h_ice_f - h_ice_i
    
        #set surface temp to freezing
        temp_sw[0] = temp_swfz
        
        #find mine heat flux through ice
        h_ice_mean = h_ice_f.mean()
        F_i = -k_ice*(temp_ice_surf_i-temp_swfz)/h_ice_mean
        
        print "F_i: %.2f, F_ocean: %.2f, h_i=%.2f, dh: %.4f" %(F_i, F_ocean, h_ice_i, dh_ice)

    else:
        
        #use thin ice algorithm
        h_ice_f, temp_ice_surf_f, temp_sw, sal_sw = modify_thin_ice(h_ice_i, temp_ice_surf_i, temp_sw, sal_sw, rho_sw, F_atm, F_ocean, alpha, params)
        dh_ice = h_ice_f - h_ice_i
        # #use F_atm and F_ocean to melt ice
        # assert F_atm<=0, "Flux inconsistency. Atmosphere warmer than ice but net surface flux is negative."
        #
        # if F_atm >=0:
        #     print "melting sea ice..."
        #     dh_ice = (F_atm - F_ocean)/(L_ice*rho_ice)
        #
        # print "F_atm: %.2f, F_ocean: %.2f, h_i=%.2f, dh: %.4f" %(F_atm, F_ocean, h_ice_i, dh_ice)
        #
      
    temp_ice_surf_f = temp_ice_surf_i
    #compute salinity change  
    dsal_sw = alpha*dh_ice*(sal_sw[:bdry_lyr].mean()-sal_ice)/(dz*bdry_lyr)
    sal_sw[:bdry_lyr] = sal_sw[:bdry_lyr]+dsal_sw


    return h_ice_f, temp_ice_surf_f, temp_sw, sal_sw, F_i
            
            
        
    
    
    # if F_i - F_ocean :
    # dh_ice = (F_i - F_ocean)/(L_ice*rho_ice)
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    