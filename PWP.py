"""
This is a Python implementation of the Price Weller Pinkel (PWP) ocean mixed layer model.
This code is based on the MATLAB implementation of the PWP model originally
written by Peter Lazarevich and Scott Stoermer (http://www.po.gso.uri.edu/rafos/research/pwp/)
(U. Rhode Island) and later modified by Byron Kilbourne (University of Washington)
and Sarah Dewey (University of Washington).

The run() function provides an outline of how the code works.

Earle Wilson
School of Oceanography
University of Washington
April 18, 2016
"""

import numpy as np
import seawater as sw
import matplotlib.pyplot as plt
from IPython.core.debugger import Tracer
import xarray as xr
import pickle
import timeit
import os
from datetime import datetime
import PWP_helper as phf
import PWP_ice
import imp

imp.reload(phf)
imp.reload(PWP_ice)

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
    
    #assert np.sum(absrb)==1, "absorption profile doesn't sum to exactly 1.0"
    #Potential issue: absrb doesn't sum to exactly 1.0 (this process leaks heat)
    
    return absrb

def pwpgo(forcing, params, pwp_out, makeLivePlots=False):
    
    """
    This is the main driver of the PWP module.
    """
    print_ice_warning = True
    print_lead_warning = True
    transition_ice_frac = False
    
    #define reference salinity and density (use these when computing heat fluxes)
    sal_ref = 34.0
    dens_ref = 1026.0
    
    #unpack some of the variables (I could probably do this more elegantly)
    F_in = forcing['F_in']
    F_out = forcing['F_out']
    emp = forcing['emp']
    taux = forcing['tx']
    tauy = forcing['ty']
    absrb = forcing['absrb']
    forcing['icec2'] = forcing['icec'].copy()
    
    
    z = pwp_out['z']
    dz = pwp_out['dz']
    dt = pwp_out['dt']
    zlen = len(z)
    tlen = len(pwp_out['time'])
    
    rb = params['rb']
    rg = params['rg']
    f = params['f']
    cpw = params['cpw']
    g = params['g']
    ucon = params['ucon']
    # alpha = params['alpha']
    
    F_net = F_in-F_out
    F_net = F_net
    
    #add dz and dt to params
    params['dt'] = dt
    params['dz'] = dz
    
    #create output variable for ocean and atmospheric heat flux 
    #TODO: use consistent naming convention for heat fluxes - pick either F or q
    pwp_out['F_oi'] = np.zeros(len(pwp_out['ice_thickness'])) #ocean-ice heat flux
    pwp_out['F_ao'] = np.zeros(len(pwp_out['ice_thickness'])) #atmosphere-ocean heat flux
    pwp_out['F_ai'] = np.zeros(len(pwp_out['ice_thickness'])) #atmosphere-ice heat flux
    pwp_out['F_i'] = np.zeros(len(pwp_out['ice_thickness'])) #heat flux through the ice
    pwp_out['F_ent'] = np.zeros(len(pwp_out['ice_thickness']))
    
    pwp_out['F_sens_ao'] = np.ma.masked_all(len(pwp_out['ice_thickness']))
    pwp_out['F_lat_ao'] = np.ma.masked_all(len(pwp_out['ice_thickness']))
    pwp_out['F_lw_ao'] = np.ma.masked_all(len(pwp_out['ice_thickness']))
    pwp_out['F_net_ao'] = np.ma.masked_all(len(pwp_out['ice_thickness']))
    
    pwp_out['F_sens_ai'] = np.ma.masked_all(len(pwp_out['ice_thickness']))
    pwp_out['F_lat_ai'] = np.ma.masked_all(len(pwp_out['ice_thickness']))
    pwp_out['F_lw_ai'] = np.ma.masked_all(len(pwp_out['ice_thickness']))
    pwp_out['F_net_ai'] = np.ma.masked_all(len(pwp_out['ice_thickness']))
    
    pwp_out['mld_exact'] = np.ma.masked_all(len(pwp_out['mld']))
    pwp_out['mld_exact2'] = np.ma.masked_all(len(pwp_out['mld']))
    pwp_out['mlt'] = np.ma.masked_all(len(pwp_out['mld']))
    pwp_out['mls'] = np.ma.masked_all(len(pwp_out['mld']))
    pwp_out['mlt_elev'] = np.ma.masked_all(len(pwp_out['mld']))
    pwp_out['alpha_true'] = np.ma.masked_all(len(pwp_out['mld']))
    pwp_out['ice_start'] = []
    
    #initialize ice thickness:
    pwp_out['ice_thickness'][0] = params['h_i0']
    if params['h_i0'] > 0.0:
        print("Initializing ice model with %sm slab of ice." %params['h_i0'])

    
    print("Number of time steps: %s" %tlen)
    
    if makeLivePlots:
        #make plot of initial profile
        phf.livePlots(pwp_out, 0)
    
    
    for n in range(1,tlen):

        print('==============================================')
        percent_comp = 100*n/float(tlen)
        yr = np.floor(pwp_out['time'][n]/365)+1; day =pwp_out['time'][n]%365
        print('Loop iter. %s (%.1f %%). Year: %i, Day: %.2f' %(n, percent_comp, yr, day))
        
        #note: n is one time step into the future, n-1 is the present

        
        #select profile data from last time step
        temp = pwp_out['temp'][:, n-1].copy()
        sal = pwp_out['sal'][:, n-1].copy()
        dens = pwp_out['dens'][:, n-1].copy()
        uvel = pwp_out['uvel'][:, n-1].copy()
        vvel = pwp_out['vvel'][:, n-1].copy()
        ps = pwp_out['ps'][:, n-1].copy()
        
        h_ice = pwp_out['ice_thickness'][n-1].copy()
        temp_ice_surf = pwp_out['surf_ice_temp'][n-1].copy()
        
        
        #compute mlt and mls before applying surface fluxes
        mld_pre, mld_idx_pre = getMLD(dens, z, params)
        
        
        #select previous forcing data (TODO: implement switch for ice conc.)
        if params['fix_alpha']:
            alpha_n = params['alpha']
        else:
            alpha_n = forcing['icec2'][n-1]
        
        skt_n = forcing['skt'][n-1]
        
        
        ### Absorb solar radiation and FWF in surf layer ###
        if h_ice==0:
            
            #set alpha (ice fraction) to 0.0:
            alpha_n = 0.0
            
            if params['use_Bulk_Formula'] == True:
                #computes atmosphere-ocean fluxes - everything but shortwave
                F_lw_ao, F_sens_ao, F_lat_ao = get_atm_ocean_HF(temp[0], forcing, alpha_n, n)
                F_sens_ao = F_sens_ao + params['qnet_offset'] #add artificial flux (default zero) to adjust sensible flux. Use to speed up cool.
                F_in_ao = (1-alpha_n)*F_in[n-1]
                F_out_ao = -(F_lw_ao+F_sens_ao+F_lat_ao) #fluxes were defined as positive down. ice fraction already accoutned for
                F_net_ao = F_in_ao-F_out_ao
                
                pwp_out['F_sens_ao'][n-1] = F_sens_ao
                pwp_out['F_lat_ao'][n-1] = F_lat_ao
                pwp_out['F_lw_ao'][n-1] = F_lw_ao
                pwp_out['F_net_ao'][n-1] = F_net_ao
            
            else:
                F_in_ao = F_in[n-1]
                F_out_ao = F_out[n-1]+ params['qnet_offset']
                F_net_ao = F_net[n-1]
            
            
            #update layer 1 temp and sal
            temp[0] = temp[0] + (F_in_ao*absrb[0]-F_out_ao)*dt/(dz*dens_ref*cpw)
            sal[0] = sal[0] + sal_ref*emp[n-1]*dt/dz
            
            ### Absorb rad. at depth ###
            temp[1:] = temp[1:] + F_in_ao*absrb[1:]*dt/(dz*dens_ref*cpw)
            
            #check if temp is less than freezing point
            T_fz = sw.fp(sal_ref, p=dz) 
            ice_heating = 0.0
            if temp[0] < T_fz:

                ice_heating = (T_fz-temp[0])*dens_ref*cpw*dz/dt #need to add artificial warming flux to compensate for ice growth
                if params['ice_ON']:
                    pwp_out['ice_start'].append(pwp_out['time'][n])
                    #generate sea ice (aka frazil ice)
                    h_ice, temp_ice_surf, temp, sal = PWP_ice.create_initial_ice(h_ice, temp_ice_surf, temp, sal, dens, alpha_n, params)
                    pwp_out['surf_ice_temp'][n] = temp_ice_surf
                    pwp_out['ice_thickness'][n] = h_ice

                else:
                    temp[0] = T_fz
                    if print_ice_warning:
                        print("surface has reached freezing temp. However, ice creation is either turned off or ice_conc is set to zero.")
                        print_ice_warning = False
            
            pwp_out['F_ao'][n-1] = F_net_ao+ice_heating
        
        else:
            
            
            #todo: implement a smooth transition in ice-fraction from open water to non-zero ice percentage
            #without this, ice frac can abruptly transition from open ocean to >50% ice cover
            
            # if transition_ice_frac:
            #     alpha_true = pwp_out['alpha_true'][n-1]
            #     alpha_forc = alpha_n
            #     d_alpha = alpha_true-alpha_forc
            #     if np.abs(d_alpha)>0.05:
            #         t_adj = int(10*(1./params['dt_d'])) #give model 10 days to catch up to the ice frac forcing
            #         alpha_nt = forcing['icec2'][n-2+t_adj]
            #         alpha_adj = np.linspace(alpha_true, alpha_nt, t_adj)
            #         forcing['icec2'][n-2:n-2+t_adj] = pwp_out['alpha_true'][n-2:n-2+t_adj]+alpha_adj
            #         alpha_n = forcing['icec2'][n-1]
            #         transition_ice_frac = False
            #
            #         #debug_here()
            
            if params['ice_ON']:
                
                if params['use_Bulk_Formula'] == True:
                    
                    #computes everything but shortwave
                    F_lw_ao, F_sens_ao, F_lat_ao = get_atm_ocean_HF(temp[0], forcing, alpha_n, n)
                    F_in_ao = (1-alpha_n)*F_in[n-1]
                    F_out_ao = -(F_lw_ao+F_sens_ao+F_lat_ao) #fluxes were defined as positive down. ice fraction already accoutned for
                    F_net_ao = F_in_ao-F_out_ao
                    
                    pwp_out['F_sens_ao'][n-1] = F_sens_ao
                    pwp_out['F_lat_ao'][n-1] = F_lat_ao
                    pwp_out['F_lw_ao'][n-1] = F_lw_ao
                    pwp_out['F_net_ao'][n-1] = F_net_ao
                    
                    F_lw_ai, F_sens_ai, F_lat_ai = get_atm_ice_HF(skt_n, forcing, alpha_n, n)
                    F_in_ai = alpha_n*F_in[n-1]
                    F_out_ai = -(F_lw_ai+F_sens_ai+F_lat_ai) #fluxes were defined as positive down. ice fraction already accoutned for
                    F_net_ai = F_in_ai - F_out_ai
                    F_ai = F_net_ai
                    
                    pwp_out['F_sens_ai'][n-1] = F_sens_ai
                    pwp_out['F_lat_ai'][n-1] = F_lat_ai
                    pwp_out['F_lw_ai'][n-1] = F_lw_ai
                    pwp_out['F_net_ai'][n-1] = F_net_ai
                
                else:
                    F_in_ao = (1-alpha_n)*F_in[n-1]
                    F_out_ao = (1-alpha_n)*F_out[n-1]
                    F_net_ao = F_in_ao - F_out_ao
                    
                    F_in_ai = alpha_n*F_in[n-1]
                    F_out_ai = alpha_n*F_out[n-1]
                    F_net_ai = F_in_ai - F_out_ai
                    F_ai = F_net_ai
                
                print("F_ao: %.2f" %F_net_ao)
                
                
                # F_out_ao = F_out_ao + (1-alpha_n)*params['qnet_offset']
                # F_out_ai = F_out_ai + alpha_n*params['qnet_offset']
                
                #compute ocean->ice heat flux
                F_oi = PWP_ice.get_ocean_ice_heat_flux(temp, sal, dens, params)
                
                #modify existing sea ice
                h_ice, temp_ice_surf, temp, sal, F_i, F_aio = PWP_ice.ice_model_T(h_ice, skt_n, temp, sal, dens, F_ai, F_oi, alpha_n, params)
                
                print(temp[:15].mean())
                
                # apply E-P flux through leads
                sal[0] = sal[0] + sal_ref*(1-alpha_n)*emp[n-1]*dt/dz #TODO: keep track of rain/snow on top of ice (big task)
                
                # apply heat flux through leads
                temp = temp + F_in_ao*absrb*dt/(dz*dens_ref*cpw) #incoming solar
                temp[0] = temp[0] - F_out_ao*dt/(dz*dens_ref*cpw) #outgoing heat
                
                #TODO: apply passive scalar flux through leads
                
                #check if temp is less than freezing point
                T_fz = sw.fp(sal_ref, p=dz)  #is sal[0] appropriate here? This is essentially brine water
                dT = temp[0]-T_fz
                lead_ice_heating = 0.0
                if dT<0:
                    print("Creating frazil ice in leads.")
                    lead_ice_heating = -dT*dens_ref*cpw*dz/dt #artificial warming flux to compensate for ice growth
                    #generate sea ice (aka frazil ice)
                    h_ice_lead, temp_ice_surf_lead, temp, sal = PWP_ice.create_initial_ice(0.0, np.nan, temp, sal, dens, (1-alpha_n), params)
                    # pwp_out['surf_ice_temp'][n] = temp_ice_surf
                    h_ice = h_ice+h_ice_lead*(1-alpha_n) #add lead ice on top of existing ice (to conserve salt/FW)
                    
                
                #save ice related output
                pwp_out['surf_ice_temp'][n] = temp_ice_surf
                pwp_out['ice_thickness'][n] = h_ice
                pwp_out['F_oi'][n-1] = F_oi
                pwp_out['F_ai'][n-1] = F_ai-F_aio
                pwp_out['F_ao'][n-1] = F_net_ao+F_aio+lead_ice_heating
                pwp_out['F_i'][n-1] = F_i
                pwp_out['alpha_true'][n] = alpha_n
            
            else:
                
                print("Whoops! This is unexpected. ")
                print("Need to turn on ice physics.")
                debug_here()
        
        #make sure ocean temp change is consistent with applied heating
        # col_mean_dQ = np.mean(temp-pwp_out['temp'][:, n-1])*dens_ref*cpw*(z.max()+dz)/dt
        # F_net = pwp_out['F_ao'][n-1]-pwp_out['F_oi'][n-1]
        # dT_error = np.abs(col_mean_dQ-F_net)
        # if dT_error>0.01:
        #     print("Warning: heat applied is not exactly consistent with heat change in water column")
        #     #debug_here()
    
            
        pwp_out['alpha_true'][n-1] = alpha_n
        
        ### compute new density ###    
        dens = getDensity(sal, temp, z, dopt=params['dens_option'])
        
        
        #save pre-entrainment temp (just for debugging)
        temp_pre = temp.copy()
        sal_pre = sal.copy()
        dens_pre = dens.copy()
        
        ### relieve static instability ###
        temp, sal, dens, uvel, vvel, ps = remove_si(z, temp, sal, dens, uvel, vvel, ps, params)
        
        if np.any(np.diff(dens)<0):
            debug_here()
        
        #compute mlt and mls after entrainment (just for debugging)
        mld_post_ent, mld_idx_post_ent = getMLD(dens, z, params)
        dMLD = mld_post_ent - mld_pre
        
        if dMLD<=0:
            F_ent = 0
        else:
            #compute change in heat content of entrained layer
            dT = np.abs(temp[mld_idx_pre:mld_idx_post_ent].mean() - temp_pre[mld_idx_pre:mld_idx_post_ent].mean())
            F_ent = dMLD*dens_ref*cpw*dT/params['dt']
        
        ### Compute MLD ###
        mld, mld_idx = getMLD(dens, z, params)
        
        ### Rotate u,v do wind input, rotate again, apply mixing ###
        ang = -f*dt/2
        uvel, vvel = rot(uvel, vvel, ang)
        du = (taux[n-1]/(mld*dens[0]))*dt
        dv = (tauy[n-1]/(mld*dens[0]))*dt
        uvel[:mld_idx] = uvel[:mld_idx]+du
        vvel[:mld_idx] = vvel[:mld_idx]+dv

        ### Apply drag to current ###
        #Original comment: this is a horrible parameterization of inertial-internal wave dispersion
        #I don't know why this is horrible -EW
        if params['drag_ON']:
            if ucon > 1e-10:
                uvel = uvel*(1-dt*ucon)
                vvel = vvel*(1-dt*ucon)
        
        uvel, vvel = rot(uvel, vvel, ang)
        
        dens1 = dens.copy()
        if np.any(np.diff(dens1)<0):
            debug_here()
        
        ### Apply Bulk Richardson number instability form of mixing (as in PWP) ###
        if rb > 1e-5:
            temp, sal, dens, uvel, vvel, ps = bulk_mix(z, temp, sal, dens, uvel, vvel, ps, dz, g, rb, zlen, mld_idx, params)
         
        dens2 = dens.copy()
        if np.any(np.diff(dens2)<0):
            debug_here()  
        
        ### Do the gradient Richardson number instability form of mixing ###
        if rg > 0:
            #TODO: get rid off extra inputs that can be passed through params.
            temp, sal, dens, uvel, vvel, ps = grad_mix(z, temp, sal, dens, uvel, vvel, ps, dz, g, rg, zlen, params, n)
            
        dens3 = dens.copy()
        if np.any(np.diff(dens3)<0):
            debug_here()
            
        ### Apply diffusion ###
        if params['rkz'] > 0:
            
            ##################################################
            """
            Here I am experimenting with applying diffusion to just the upper portion of the water column.
            By default, diffusion is applied everywhere. For long (multi-year) model runs, this can lead to significant 
            smoothing in the deep ocean. The params['diff_zlim'] sets the maximum depth over which
            diffusion is applied - this max depth must be below the MLD. Exercisiing this option preserves 
            the deep ocean stratification but it comes at a cost of higher numerical error.
            """
            nz1 = len(z[z<=params['diff_zlim']]) #by default is params['diff_zlim'] is some big number
            
            #Limit diffusion to some upper level of ocean but no shallower than the maximum MLD
            nz2 = mld_idx
            if nz2>nz1:
                print("WARNING: MLD (%.1fm) exceeds diffusion z-limit (%.1fm)." %(mld, params['diff_zlim']))
                print("Resetting diffusion z-limit to %.1fm." %(z[nz2+2]))
                nz=nz2+2
            else:
                nz=nz1
    
            ##################################################
            
            
            temp = diffus(params['dstab'], nz, temp)
            sal = diffus(params['dstab'], nz, sal)
            uvel = diffus(params['dstab'], nz, uvel)
            vvel = diffus(params['dstab'], nz, vvel)
            ps = diffus(params['dstab'], nz, ps)
            
            ### compute new density ###    
            dens = getDensity(sal, temp, z, dopt=params['dens_option'])
                
            

        
        # find MLD again, after all mixing processes complete
        mld_idx = np.flatnonzero(dens-dens[0]>params['mld_thresh'])[0] #finds the first index that exceeds ML threshold
        mld_post_mix = z[mld_idx]

        # find MLD by interpolating to the exact value
        mld_idx2 = np.flatnonzero(dens-dens[0]==0)[-1] #finds the last index of ML
        from scipy.interpolate import interp1d
        p_int = interp1d(dens[mld_idx2:], z[mld_idx2:]) # interp1d does not work well with perfectly homogenous mixed layers
        mld_exact = p_int(dens[mld_idx2]+params['mld_thresh'])
        mld_exact2 = p_int(dens[mld_idx2]+0.01) #
        
        #save MLT, MLS and MLT_elevation
        MLS_post = np.mean(sal[:mld_idx]);
        MLT_post = np.mean(temp[:mld_idx])
        #T_fz = sw.fp(MLS_post, 0)
        T_elev = MLT_post - T_fz
    
        ### update output profile data ###
        pwp_out['temp'][:, n] = temp
        pwp_out['sal'][:, n] = sal
        pwp_out['dens'][:, n] = dens
        pwp_out['uvel'][:, n] = uvel
        pwp_out['vvel'][:, n] = vvel
        pwp_out['ps'][:, n] = ps
        pwp_out['mld'][n] = mld_post_mix
        pwp_out['mlt'][n] = MLT_post
        pwp_out['mls'][n] = MLS_post
        pwp_out['mlt_elev'][n] = T_elev
        pwp_out['mld_exact'][n] = mld_exact
        pwp_out['mld_exact2'][n] = mld_exact2
        # pwp_out['F_atm'][n] = (1-alpha_n)*F_atm
        pwp_out['F_ent'][n] = F_ent
        pwp_out['alpha_true'][n] = alpha_n
        
        #show makeLivePlots
        #debug_here()
        if makeLivePlots:
            phf.livePlots(pwp_out, n)
        
    pwp_out['ice_start'] = np.array(pwp_out['ice_start'])
    
    return pwp_out    

def remove_si(z, t, s, d, u, v, ps, params):
    
    # Find and relieve static instability that may occur in the
    # density array 'd'. This simulates free convection.
    # ml_index is the index of the depth of the surface mixed layer after adjustment,
    
    stat_unstable = True
    
    while stat_unstable:
        
        d_diff = np.diff(d)
        if np.any(d_diff<0): #maybe use a small tolerance??
            stat_unstable=True
            first_inst_idx = np.flatnonzero(d_diff<0)[0]+1 #+1 because of diff function
            #d0 = d
            (t, s, d, u, v, ps) = mix5(z, t, s, d, u, v, ps, first_inst_idx, params)
            
            #debug_here()
            
            #plot density
            # plt.figure(num=86)
            # plt.clf() #probably redundant
            # plt.plot(d0-1000, range(len(d0)), 'b-', label='pre-mix')
            # plt.plot(d-1000, range(len(d0)), 'r-', label='post-mix')
            # plt.gca().invert_yaxis()
            # plt.xlabel('Density-1000 (kg/m3)')
            # plt.grid(True)
            # plt.pause(0.05)
            # plt.show()
            
        else:
            stat_unstable = False
    
    return t, s, d, u, v, ps    

def mix5(z, t, s, d, u, v, ps, j, params):
    
    #This subroutine mixes the arrays t, s, u, v down to level j (level where there is instability).
    j = j+1 #so that the j-th layer is included in the mixing
    t[:j] = np.mean(t[:j])
    s[:j] = np.mean(s[:j])
    
    ### compute new density ###    
    d[:j]  = getDensity(s[:j], t[:j], z[:j], dopt=params['dens_option'])

    #d[:j] = np.mean(d[:j]) #doing ignores the non-linearities in the eqn of state
    u[:j] = np.mean(u[:j])
    v[:j] = np.mean(v[:j])
    ps[:j] = np.mean(ps[:j])
    
    return t, s, d, u, v, ps

def rot(u, v, ang):
    
    #This subroutine rotates the vector (u,v) through an angle, ang
    r = (u+1j*v)*np.exp(1j*ang)
    u = r.real
    v = r.imag
    
    return u, v

def bulk_mix(z, t, s, d, u, v, ps, dz, g, rb, nz, ml_idx, params):
    #sub-routine to do bulk richardson mixing
    
    rvc = rb #critical rich number??
    
    for j in range(ml_idx, nz):
        h   = z[j]
        dd  = (d[j]-d[0])/d[0]
        dv  = (u[j]-u[0])**2+(v[j]-v[0])**2
        
        if dv == 0:
            rv = np.inf
        else:
            rv = g*h*dd/dv
        
        if rv > rvc:
            break
        
        else:
            t, s, d, u, v, ps = mix5(z, t, s, d, u, v, ps, j, params)
    
    return t, s, d, u, v, ps

def grad_mix(z, t, s, d, u, v, ps, dz, g, rg, nz, params, n):
    
    #TODO: get rid off function arguments that are already in params
    
    #copied from source script:
    # %  This function performs the gradeint Richardson Number relaxation
    # %  by mixing adjacent cells just enough to bring them to a new
    # %  Richardson Number.
    # %  Compute the gradeint Richardson Number, taking care to avoid dividing by
    # %  zero in the mixed layer.  The numerical values of the minimum allowable
    # %  density and velocity differences are entirely arbitrary, and should not
    # %  effect the calculations (except that on some occasions they evidently have!)
    #print "entered grad mix"
    
    rc = params['rg'] #critical rich. number
    j1 = 0
    j2 = nz-1
    j_range = np.arange(j1, j2)
    i = 0 #loop count
    #debug_here()
    
    while 1:
        #TODO: find a better way to do implement this loop
        
        r = np.zeros(len(j_range),)
        
        for j in j_range:
            dd = (d[j+1]-d[j])/d[j]
            dv = (u[j+1]-u[j])**2+(v[j+1]-v[j])**2
            if dv == 0:
                r[j] = np.inf
            else:
                #compute grad. rich. number
                r[j] = params['g']*params['dz']*dd/dv
            
            # if r[j]<0:
            #     print("whoa there. negative richardson number")
            #     print(r[j])
            #     debug_here()
        
        #find the smallest value of r in the profile
        r_min = np.min(r)
        j_min_idx = np.argmin(r)
        
        #Check to see whether the smallest r is above critical or not.
        if r_min > rc:
            break
        
        #Mix the cells j_min_idx and j_min_idx+1 that had the smallest Richardson Number
        t, s, d, u, v, ps = stir(z, t, s, d, u, v, ps, rc, r_min, j_min_idx, params, n)
        
        #recompute the rich number over the part of the profile that has changed
        j1 = j_min_idx-2
        if j1 < 1:
             j1 = 0
        
        j2 = j_min_idx+2
        if j2 > nz-1:
             j2 = nz-1
        
        i+=1
        
        #Check to make sure profile contains no nans
        if np.any(np.isnan(t)) or np.any(np.isnan(s)) or np.any(np.isnan(d)):
            debug_here()
            raise ValueError("Nans appeared in profile. Something very bad happened!")
        
    
    return t, s, d, u, v, ps

def stir(z, t, s, d, u, v, ps, rc, r, j, params, n):
    
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
    rcon = 0.02+(rc-r)/2.
    rnew = rc+rcon/5.
    f = 1-r/rnew
    
    #mix temp
    dt = (t[j+1]-t[j])*f/2.
    t[j+1] = t[j+1]-dt
    t[j] = t[j]+dt
    
    #mix sal
    ds = (s[j+1]-s[j])*f/2.
    s[j+1] = s[j+1]-ds
    s[j] = s[j]+ds
    
    #mix passive scalar
    d_ps = (ps[j+1]-ps[j])*f/2.
    ps[j+1] = ps[j+1]-d_ps
    ps[j] = ps[j]+d_ps
    
    #recompute density
    #have to be careful here. x[j:j+1] in python is not the same as x[[j,j+1]]. We want the latter
    d[[j,j+1]] = getDensity(s[[j,j+1]], t[[j,j+1]], z[[j,j+1]], dopt=params['dens_option'])

    # if n==512:
    #     debug_here()
    
    #mix velocities
    du = (u[j+1]-u[j])*f/2.
    u[j+1] = u[j+1]-du
    u[j] = u[j]+du
    
    dv = (v[j+1]-v[j])*f/2.
    v[j+1] = v[j+1]-dv
    v[j] = v[j]+dv
    
    return t, s, d, u, v, ps

def diffus(dstab,nz,a):
    
    "finite difference implementation of diffusion equation"
    
    #matlab code:
    #a(2:nz-1) = a(2:nz-1) + dstab*(a(1:nz-2) - 2*a(2:nz-1) + a(3:nz));
    
    a[1:nz-1] = a[1:nz-1] + dstab*(a[0:nz-2] - 2*a[1:nz-1] + a[2:nz])
    return a

def getMLD(dens, z, params):
    
    #find ml index
    mld_idx = np.flatnonzero(dens-dens[0]>params['mld_thresh'])[0]
    
    #check to ensure that ML is defined
    assert mld_idx.size is not 0, "Error: Mixed layer depth is undefined."
    
    #get MLD
    mld = z[mld_idx-1]
    
    if mld_idx==1:
        mld = params['dz']
    
    return mld, mld_idx

def get_atm_ocean_HF(sst, forcing, alpha, n, height_match=False):
    
    """
    Function to compute ocean atmosphere heat fluxes in the presence of sea ice.
    
    flux sign convention: postive-down
    
    This is based on Large and Yeager 2004
    Q = f_i*Q_io + f_o*Q_ao
    """
    
    def psi(zeta):
        
        X_o = (1-16*zeta)**(0.25)
        
        if zeta>=0:
            #stable
            psi_m = -5*zeta
            psi_h = -5*zeta
        else:
            #unstable
            psi_m = 2*np.log((1.+X_o)/2.) + np.log((1.+X_o**2)/2.) - 2*np.arctan(X_o) + np.pi/2
            psi_h = 2*np.log((1.+X_o**2)/2.)
        
        return psi_m, psi_h
    
    
    sst_K = sst+273.13
    
    #define params
    f_i = alpha
    f_o = 1-f_i
    q1_ocn = 0.98*640380. #kg/m3 (Large and Yeager 2004 pg. 7)
    q2_ocn = -5107.4 #K
    rho_air = 1.22 #kg/m3
    c_p_air = 1000.5 #J/kg/K (specific heat of air)
    L_v = 2.5e6 #J/kg (Latent heat of evaporation)
    kappa = 0.4 #von karman constant
    g = 9.81 #acc due to gracity
    r = 287 #J/kg/K ideal gas constant
    
    #get key atmospheric forcing variables
    q_2m = forcing['shum2m'][n-1]
    atemp2m = forcing['atemp2m'][n-1]
    u10m = forcing['u10m'][n-1]
    v10m = forcing['v10m'][n-1]
    tx = forcing['tx'][n-1]
    ty = forcing['ty'][n-1]
    
    #get longwave over the ocean
    sigma = 5.67e-8 #W/m2/K4 (step-boltzmann constant)
    eps = 1.0 #emissivity
    F_lw_ao = f_o*(forcing['dlw'][n-1] - eps*sigma*sst_K**4)
    
    #get variables for turbulent heat fluxes
    Un_10m_init = np.sqrt(u10m**2 + v10m**2) #10m wind speed
    q_sat_ocn = (1./rho_air)*q1_ocn*np.e**(q2_ocn/sst_K) #saturation specific humidity
    theta_2m = atemp2m*(1000./1002.)**(r/c_p_air) #potential air temp
    theta_v_2m = theta_2m*(1. + 0.608*q_2m) #virtual potential air temp
    
    #define turbulent heat transfer coefficients (assuming unstable atmosphere)
    C_d = (1e-3)*(2.7/Un_10m_init + 0.142 + Un_10m_init/13.09)
    C_e = (1e-3)*34.6*np.sqrt(C_d)
    C_h = (1e-3)*32.7*np.sqrt(C_d)
    
    # the next series of steps are from Large and Yeager (2004) pg.8-9
    # the goal is to shift the air temp and spec. humidity from 2m to 10m
    
    #define a few more parameters
    # u_star = np.sqrt(np.sqrt(tx**2 + ty**2)/rho_air)
    u_star = np.sqrt(C_d)*Un_10m_init
    # theta_star_2m_ocn = theta_2m - sst_K
    # q_star_2m_ocn = q_2m - q_sat_ocn
    t_star_2m_ocn = (C_h/np.sqrt(C_d))*(theta_2m - sst_K)
    q_star_2m_ocn = (C_e/np.sqrt(C_d))*(q_2m - q_sat_ocn)
    
    
    #compute initial sensible and latent heat fluxes
    F_sens_ao_init = f_o*rho_air*c_p_air*C_h*(theta_2m - sst_K)*Un_10m_init
    F_lat_ao_init = f_o*L_v*rho_air*C_e*(q_2m - q_sat_ocn)*Un_10m_init
    
    if height_match==False:
        
        if n==1:
            print("Warning: Skipping height match between winds and humidity!")
            
        print("F_sens_ao = %.2f, F_lat_ao = %.2f, F_lw_ao = %.2f W/m2" %(F_sens_ao_init, F_lat_ao_init, F_lw_ao))
        
        return F_lw_ao, F_sens_ao_init, F_lat_ao_init
        
    #WARNING: The code below sometimes produces really weird values. Generally, the adjustments are relatively small, tending to slightly 
    #increase the heat fluxes. However, on rare occassions these differences can be very large O(1000 W/m2). This is obviously unphysical. It's probably best
    #to ignore the height adjustment process for now and stick to the simple bulk forumla computations.
    
        
    #get monin-obukhov length scale
    # L_mo_ocn = -rho_air*c_p_air*theta_v_2m*u_star**3/(kappa*g*qsens_ao_init)
    
    #compute stability parameters
    z_t = 2.
    z_q = 2.
    z_u = 10.
    zeta_2m = (kappa*g*z_t/u_star**2)*(t_star_2m_ocn/theta_v_2m + q_star_2m_ocn/(q_2m + 1./0.608))
    zeta_10m = (kappa*g*z_u/u_star**2)*(t_star_2m_ocn/theta_v_2m + q_star_2m_ocn/(q_2m + 1./0.608))
    
    psi_m_2m, psi_h_2m  = psi(zeta_2m)
    psi_m_10m, psi_h_10m  = psi(zeta_10m)
    
    #shift wind speed to 10m and neutral stability, and temp and hum to the wind height
    Un_10m = Un_10m_init*(1 + (np.sqrt(C_d)/kappa)*(np.log(z_u/z_u) - psi_m_10m))**(-1)
    theta_10m = theta_2m - (t_star_2m_ocn/kappa)*(np.log(z_t/z_u) + psi_h_10m - psi_h_2m) #paper says t_star instead of q_star???
    q_10m = q_2m - (q_star_2m_ocn/kappa)*(np.log(z_t/z_u) + psi_h_10m - psi_h_2m)
    
    #update intial transfer coefficients
    C_d_new = C_d*(1+ (np.sqrt(C_d)/kappa)*(np.log(z_u/z_u) - psi_m_10m) )**(-2.)
    C_h_new = C_h*np.sqrt(C_d_new/C_d)*(1+ (C_h/(kappa*np.sqrt(C_d)) )*(np.log(z_u/z_u) - psi_h_10m))**(-1.)
    C_e_new = C_e*np.sqrt(C_d_new/C_d)*(1+ (C_e/(kappa*np.sqrt(C_d)) )*(np.log(z_u/z_u) - psi_h_10m))**(-1.)
    
    #recompute qsens and qlat
    # theta_star_10m_ocn = theta_10m - sst_K
    # q_star_10m_ocn = q_10m - q_sat_ocn
    
    F_sens_ao = f_o*rho_air*c_p_air*C_h_new*(theta_10m - sst_K)*Un_10m #if atm warmer than ocean F_sens is positive (warms ocean)
    F_lat_ao = f_o*L_v*rho_air*C_e_new*(q_10m - q_sat_ocn)*Un_10m #if atm warmer than ocean F_sens is positive (warms ocean)
    

    
    # theta_v_10m = theta_10m*(1. + 0.608*q_10m)
    # t_star_10m_ocn = (C_h/np.sqrt(C_d))*theta_star_10m_ocn
    # q_star_10m_ocn = (C_e/np.sqrt(C_d))*q_star_10m_ocn
    
    if abs(F_sens_ao)>1000:
        print("woah!")
        debug_here()
    
    
    print("F_sens_ao = %.2f, F_lat_ao = %.2f, F_lw_ao = %.2f W/m2" %(F_sens_ao, F_lat_ao, F_lw_ao))
    #
    # if n%500==0 and n<=1500:
    #     debug_here()
    
    return F_lw_ao, F_sens_ao, F_lat_ao
    

def get_atm_ice_HF(surf_temp, forcing, alpha, n):
    
    """
    Like get_ocean_atm_HF() but for ice-atmosphere fluxes.
    
    TODO: combine the two functions
    """
    
    surf_temp_K = surf_temp+273.13
    
    
    #define params
    f_i = alpha
    f_o = 1-f_i
    q1 = 11637800 #kg/m3 (Large and Yeager 2004 pg. 16)
    q2 = -5897.8 #K
    rho_air = 1.22 #kg/m3
    c_p_air = 1000.5 #J/kg/K (specific heat of air)
    L_v = 2.5e6 #J/kg
    L_f = 2.839e6 #J/kg (Latent heat of sublimation)
    kappa = 0.4 #von karman constant
    g = 9.81 #acc due to gracity
    r = 287 #J/kg/K ideal gas constant
    C_e,C_h,C_d = 1.63e-3, 1.63e-3, 1.63e-3
    
    #get key atmospheric forcing variables
    q_2m = forcing['shum2m'][n-1]
    atemp2m = forcing['atemp2m'][n-1]
    u10m = forcing['u10m'][n-1]
    v10m = forcing['v10m'][n-1]
    tx = forcing['tx'][n-1]
    ty = forcing['ty'][n-1]
    
    U_10m = np.sqrt(u10m**2 + v10m**2) #10m wind speed
    q_sat_ice = (1/rho_air)*q1*np.e**(q2/surf_temp_K) #saturation specific humidity
    theta_2m = atemp2m*(1000./1002.)**(r/c_p_air) #potential air temp
    theta_v_2m = theta_2m*(1. + 0.608*q_2m) #virtual potential air temp
    
    #compute latent and sens heat fluxes over ice
    F_lat_ai = f_i*L_f*rho_air*C_e*(q_2m-q_sat_ice)*U_10m
    F_sens_ai = f_i*rho_air*C_h*(theta_2m-surf_temp_K)*U_10m
    
    #compute lw radiation over ice
    sigma = 5.67e-8 #W/m2/K4 (step-boltzmann constant)
    eps = 0.95 #emissivity
    F_lw_ai = f_i*eps*(forcing['dlw'][n-1] - sigma*surf_temp_K**4)
    
    
    return F_lw_ai, F_sens_ai, F_lat_ai
    

def getDensity(s,t,z, dopt):
    
    
 
    if dopt == 'dens0':
        dens = sw.dens0(s, t)
        
    elif dopt == 'dens':
        dens = sw.dens(s, t, z)
        
    elif dopt == 'pdens':
        dens = sw.pden(s, t, z, pr=0)
    else:
        print("Error: set appropriate density option: 'dens0', 'pdens' or 'dens' ")
        debug_here()
    
    return dens
    

def local_stir(z, s, t, ps, dopt):
    
    
    s0 = s.copy()
    t0 = t.copy()
    ps0 = ps.copy()
    d = getDensity(s,t,z, dopt)
    d0 = d.copy()

    dd = np.diff(d0)
    mix_frac = 0.5
    
    print("Stabilizing profile")
    i=0
    while any(dd<0):
        
        i=i+1
        
        neg_idx = np.flatnonzero(dd<=0)
        neg_idx_r = neg_idx[::-1] #mix from bottom up
        
        for idx in neg_idx_r:
            
            j = idx
            
            #mix temp
            dt = (t[j+1]-t[j])*mix_frac/2.
            t[j+1] = t[j+1]-dt
            t[j] = t[j]+dt
    
            #mix sal
            ds = (s[j+1]-s[j])*mix_frac/2.
            s[j+1] = s[j+1]-ds
            s[j] = s[j]+ds
    
            #mix passive scalar
            d_ps = (ps[j+1]-ps[j])*mix_frac/2.
            ps[j+1] = ps[j+1]-d_ps
            ps[j] = ps[j]+d_ps
            
        
            #recompute density
            d[[j,j+1]] = getDensity(s[[j,j+1]], t[[j,j+1]], z[[j,j+1]], dopt=dopt)
            dd = np.diff(d)
        
        if mix_frac<1:
            mix_frac = mix_frac+0.1
        else:
            mix_frac=1
        
        if i>5e3:
            print("Failed to stabilize profile after %i iterations :(" %i)
            break
        
    print("Profile stabilized. %i iterations required." %i)
        
    #plot stablized profile
    fig, axes = plt.subplots(1,3, figsize=(12,6), sharey=True)
    axes[0].plot(t0, z, label='initial')
    axes[0].plot(t, z, label='stabilized')
    axes[0].grid(True)
    axes[0].invert_yaxis()
    axes[0].set_ylabel("Depth (m)")
    axes[0].set_xlabel("Temperature (C)")
    axes[0].legend(loc=0)
    axes[0].set_xlim(-2, 1.5)
    
    axes[1].plot(s0, z, label='initial')
    axes[1].plot(s, z, label='stabilized')
    axes[1].grid(True)
    axes[1].invert_yaxis()
    axes[1].set_ylabel("Depth (m)")
    axes[1].set_xlabel("Salinity (psu)")
    axes[1].legend(loc=0)
    axes[1].set_xlim(33.5, 34.75)
    
    axes[2].plot(d0-1000, z, label='initial')
    axes[2].plot(d-1000, z, label='stabilized')
    axes[2].grid(True)
    axes[2].invert_yaxis()
    axes[2].set_ylabel("Depth (m)")
    axes[2].set_xlabel("%s - 1000 (kg/m3)" %dopt)
    axes[2].legend(loc=0)
    from matplotlib.ticker import FormatStrFormatter
    axes[2].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axes[2].set_xlim(26.8, 28)
    
    debug_here()
    
    return s, t, ps, d 
        
        
    

if __name__ == "__main__":
    
    print("Running default test case using data from Beaufort gyre...")
    
    forcing_fname = 'beaufort_met.nc'
    prof_fname = 'beaufort_profile.nc'
    run(met_data=forcing_fname, prof_data=prof_fname)
