import visualization.movie_maker as mm
import numpy as np
import gadget
import cosmology as mycosmo
import scipy.misc
import scipy.optimize
import h5py
import fit_functions as fitfun
import pfh_utils as util
#import matplotlib.pylab as pylab
import sys

def test():
    #hmain = halo()
    #print hmain.a_scale
    #hmain.a_scale = 1.
    #print hmain.a_scale

    collect_main_halo_properties('../../zooms/m11_hr/tmp', restart_with_snum=400, snum_max=0 )


## class which contains all of the variables to save for the tracked halo 
##   at each timestep (add new things here, then set them in the loop routine below)
class halo:
    def __init__(self, N):
        self.fraction_of_rvir_for_inflow_shells = np.array([ 0.1, 0.25, 0.5, 1.0, 2.0 ]) # anuli for inflow
        self.temperature_threshold_list = np.array([ 1.e3, 1.e4, 1.e5, 1.e6, 1.e10 ]) # gas temperature thresholds
        self.N_profile_radii = np.array(100.)
    
        self.a_scale = np.zeros((N))
        self.time = np.zeros((N))
        self.r_200 = np.zeros((N))
        self.m_200 = np.zeros((N))
        self.v2_max_all = np.zeros((N))
        self.v2_max_dm = np.zeros((N))
        self.halo_concentration = np.zeros((N))
        self.halo_concentration_freeslope = np.zeros((N))
        self.halo_innerslope_freeslope = np.zeros((N))
        self.m_gas_halo = np.zeros((N))
        self.m_stars_halo = np.zeros((N))
        self.m_stars_new = np.zeros((N))
        self.metal_stars_halo = np.zeros((N))
        self.metal_stars_new = np.zeros((N))
        self.stellar_reff = np.zeros((N))
        self.m_baryon_new_this_timestep = np.zeros((N))
        self.metals_baryon_new_this_timestep = np.zeros((N))
        self.m_baryon_totally_new = np.zeros((N))
        self.metals_baryon_totally_new = np.zeros((N))
        
        self.center = np.zeros((N,3))
        self.com_vel = np.zeros((N,3))
        self.j_total = np.zeros((N,3))
        self.j_total_dm = np.zeros((N,3))
        self.j_stars = np.zeros((N,3))

        #self.stellar_sersic_parameters = np.zeros((N,3))
        #self.stellar_sersicdisk_parameters = np.zeros((N,6))
        
        Nshells = self.fraction_of_rvir_for_inflow_shells.size
        Ntemp = self.temperature_threshold_list.size
        
        self.mdot_inflow = np.zeros((N,Nshells))
        self.mdot_outflow = np.zeros((N,Nshells))
        self.pdot_inflow = np.zeros((N,Nshells))
        self.pdot_outflow = np.zeros((N,Nshells))
        self.mdot_inflow_metal = np.zeros((N,Nshells))
        self.mdot_outflow_metal = np.zeros((N,Nshells))
        self.mdot_inflow_hot = np.zeros((N,Nshells))
        self.mdot_outflow_hot = np.zeros((N,Nshells))
        self.pdot_inflow_hot = np.zeros((N,Nshells))
        self.pdot_outflow_hot = np.zeros((N,Nshells))
        self.mdot_inflow_metal_hot = np.zeros((N,Nshells))
        self.mdot_outflow_metal_hot = np.zeros((N,Nshells))
        
        self.m_gas_halo_tempcut = np.zeros((N,Ntemp))
        self.metal_gas_halo_tempcut = np.zeros((N,Ntemp))
        
        self.stellar_mass_profile = np.zeros((N,self.N_profile_radii))
        self.stellar_mass_profile_rmin = np.zeros((N))
        self.stellar_mass_profile_rmax = np.zeros((N))
        self.dm_mass_profile = np.zeros((N,self.N_profile_radii))
        self.dm_mass_profile_rmin = np.zeros((N))
        self.dm_mass_profile_rmax = np.zeros((N))
        self.gas_mass_profile = np.zeros((N,self.N_profile_radii))
        self.gas_mass_profile_rmin = np.zeros((N))
        self.gas_mass_profile_rmax = np.zeros((N))
        

##
## routine to loop over all the snapshots in a directory, tracking the 'main' 
##   halo and pulling various properties at each timestep, for later analysis
##
def collect_main_halo_properties( snapdir, 
    use_h0=1,four_char=0,cosmological=1,skip_bh=1, \
    restart_with_snum=0, snum_max=0):

    N_snaps = 441 # need to know this ahead of time, will look for every snap from 0-this number
    r_clip = 600. # maximum physical radius to gather initial particles
    hot_gas_virial_threshold = 0.5 # fraction of t_vir to qualify as 'hot' gas (for inflow/outflow)
    virial_threshold = 200. # standard virial overdensity threshold

    ## initialize variables needed below
    MYHALO = halo(N_snaps)
    id_baryons_in_any_halo_timestep = np.array([-1],dtype='i')
    id_baryons_in_previous_halo_timestep = np.array([-1],dtype='i')
    if(snum_max<=0): snum_max = N_snaps
    snum_min = 0

    if (restart_with_snum>0):
        ## will assume there is a file with previous snapshot info saved, and load that:
        fname = snapdir+'/halo_prop_vtime.hdf5'
        print('getting data from file : ',fname)
        infi = h5py.File(fname,'r')
        for key in infi.keys():
            MYHALO.__dict__[key] = np.copy(infi[key])
        infi.close()
        fname = snapdir+'/halo_prop_vtime_ids.hdf5'
        infi = h5py.File(fname,'r')
        id_baryons_in_previous_halo_timestep = np.copy(infi['id_baryons_in_previous_halo_timestep'])
        id_baryons_in_any_halo_timestep = np.copy(infi['id_baryons_in_any_halo_timestep'])
        infi.close()
        snum_min = restart_with_snum

    print('N-profile-radii == ',MYHALO.N_profile_radii)

    ## (for now not calling this here, but want to make sure it *has* been run before using this routine)
    ## want to define a *smooth* center, or else differential quantities can be screwed up
    #print '... building/loading camera positions list ...'
    #mm.build_camera_centering( snapshot_list, snapdir, force_rebuild=0, \
    #    center_on_bh=0,center_on_com=0,set_fixed_center=0,\
    #    use_h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh )

    print('... entering main snapshot loop ...')
    ## now enter main loop over snapshots in the simulation set
    snapshot_list = range(N_snaps+1)
    for i_snap,snum in zip(range(N_snaps),snapshot_list):
        if (snum < snum_min): continue;
        if (snum > snum_max): continue; 
        ## load header info
        print(' processing snapshot ',i_snap,' snapnumber=',snum,' in ',snapdir)  
        PPP_head = gadget.readsnap(snapdir,snum,1,header_only=1,\
            h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh)
        if (PPP_head['k']<0): continue;
        a_scale_i = PPP_head['time']
        z_redshift = 1./a_scale_i - 1.
        MYHALO.a_scale[i_snap] = a_scale_i
        MYHALO.time[i_snap] = mm.cosmological_time(a_scale_i)
        ## load center position
        cen = mm.get_precalc_zoom_center(snapdir,a_scale_i,cosmological=cosmological)
        #cen = gadget.calculate_zoom_center(snapdir,snum,cen=[0.,0.,0.],h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh)
        MYHALO.center[i_snap,:] = cen

        ## alright now load the actual --data-- to get the virial size, mass, etc.
        m,ptype,pos,vel = load_pos_vel(snapdir,snum, 
            use_h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh)
        for j in [0,1,2]: pos[:,j] -= cen[j]

        ## first step: clip to the particles of interest and             
        r2 = quadratic_sum(pos)
        ok=(r2 < r_clip*r_clip); r2=r2[ok]; m=m[ok]; ptype=ptype[ok]; r=np.sqrt(r2);
        Nok=m.size; pos_tmp=np.zeros((Nok,3)); vel_tmp=np.zeros((Nok,3));
        for j in [0,1,2]:
            pos_tmp[:,j] = pos[ok,j]
            vel_tmp[:,j] = vel[ok,j]
        ## sort them
        s=np.argsort(r); r_s=r[s]; m_s=m[s]; m_enc=np.cumsum(m_s)

        ## get the mean density inside each radius, in units of cosmic mean
        rho_background = mycosmo.mass_density( z_redshift , cgs=1)
        units_rho = 1.0e10 * 1.989e33 / ((4.*np.pi/3.)*((3.086e21)**3.)) 
        rho_mean_halo = m_enc / (r_s*r_s*r_s) * units_rho/rho_background
        
        if (rho_mean_halo[0] <= virial_threshold):
            print('No halo yet virialized detected in snap ',i_snap,j)
            continue
        
        outside = (r_s > r_clip/50.) & (rho_mean_halo < virial_threshold)
        if (r_s[outside].size <= 0):
            MYHALO.r_200[i_snap] = r_clip/50.
        else:
            MYHALO.r_200[i_snap] = r_s[outside][0]
        r_200 = MYHALO.r_200[i_snap]
            
        print('r200 = ',MYHALO.r_200[i_snap])
        
        ## clip again, now to everything inside r_vir
        ok = ( r < r_200 )
        r2=r2[ok]; r=r[ok]; m=m[ok]; ptype=ptype[ok];
        Nok=m.size; pos=np.zeros((Nok,3)); vel=np.zeros((Nok,3))
        for j in [0,1,2]:
            pos[:,j] = pos_tmp[ok,j]
            vel[:,j] = vel_tmp[ok,j]
        m_200 = m_total = np.sum(m)
        particle_wt = m / m_total
        MYHALO.m_200[i_snap] = m_total
        temp_virial = tvir_of_halo( m_200, r_200 )
        
        print('m200, log10(tvir) = ',np.log10(MYHALO.m_200[i_snap])+10., np.log10(temp_virial))
        
        ## get halo center-of-mass velocity 
        for j in [0,1,2]: MYHALO.com_vel[i_snap,j] = np.sum( particle_wt * vel[:,j] )
        ## and subtract this off so have relative velocities
        for j in [0,1,2]: vel[:,j] = vel[:,j] - MYHALO.com_vel[i_snap,j]
        
        ## get angular momentum vectors for each particle
        j_spec = np.cross( pos, vel , axis=1 )
        print('jspec : ',j_spec.shape,np.min(j_spec),np.max(j_spec),m.shape)
        for j in [0,1,2]: MYHALO.j_total[i_snap,j] = np.sum( m*j_spec[:,j] )
        print('j_total = ',MYHALO.j_total[i_snap,:])
        is_dm = (ptype != 0) & (ptype != 4)
        for j in [0,1,2]: MYHALO.j_total_dm[i_snap,j] = np.sum( m[is_dm]*j_spec[is_dm,j] )
        
        ## get vmax (for all particles and for DM profile alone)
        s=np.argsort(r); r_s=r[s]; m_s=m[s]; t_s=ptype[s]; 
        m_enc=np.cumsum(m_s); v2=m_enc/r_s;
        MYHALO.v2_max_all[i_snap] = np.max(v2)
        is_dm_s = (t_s != 0) & (t_s != 4)
        print('N tot and only dm == ',r_s.size,t_s[is_dm_s].size)
        m_enc_dm = np.cumsum(m_s[is_dm_s]); v2_dm = m_enc_dm/r_s[is_dm_s]
        MYHALO.v2_max_dm[i_snap] = np.max(v2_dm)
        
        print('v2max_all/dm = ',MYHALO.v2_max_all[i_snap],MYHALO.v2_max_dm[i_snap])
        print('vmax_all/dm = ',1.e-5*np.sqrt(MYHALO.v2_max_all[i_snap]*6.67e-8*1.989e43/3.086e21),\
            1.e-5*np.sqrt(MYHALO.v2_max_dm[i_snap]*6.67e-8*1.989e43/3.086e21))
    
        ## fit DM profile (NFW, etc)
        dlnr = 0.05*np.log(10.); rmin = r_200/30.;
        #rmin = r_200/50. ## inner slopes are quite sensitive actually to rmin, for fitting
        ## get density profile over some interpolation range 
        if (r_s[is_dm_s][0] > rmin): rmin=r_s[is_dm_s][0];
        ln_rg = np.arange( np.log(rmin), np.log(r_200), dlnr )
        rg = np.exp(ln_rg); mg = np.interp( rg, r_s[is_dm_s], m_enc_dm )
        rho_g = mg/(4.*np.pi*rg*rg*rg) * np.abs(my_simple_derivative( np.log(mg), ln_rg ))
        ln_rho_g = np.log(rho_g)
        ##pylab.plot(ln_rg/np.log(10.),ln_rho_g/np.log(10.))

        ## guess scale radius from R_e and corresponding density 
        guess_lnrscale = np.interp( m_200/2.0, mg, ln_rg )
        guess_lnpscale = np.interp( guess_lnrscale, ln_rg, ln_rho_g ) + 1.386
        ## fit NFW profile: 
        popt,pcov = scipy.optimize.curve_fit(nfw,ln_rg,ln_rho_g,p0=[guess_lnpscale,guess_lnrscale]);
        r_scale = np.exp(popt[1])
        MYHALO.halo_concentration[i_snap] = r_200 / r_scale
        ## fit NFW profile with fixed outer but free inner slope (guess NFW to start): 
        popt,pcov = scipy.optimize.curve_fit(nfw_freeinnerslope,ln_rg,ln_rho_g,\
            p0=[guess_lnpscale,guess_lnrscale,1.]);
        r_scale = np.exp(popt[1]); 
        MYHALO.halo_concentration_freeslope[i_snap] = r_200 / r_scale
        MYHALO.halo_innerslope_freeslope[i_snap] = popt[2]

        print('NFW: r/p_guess',np.exp(guess_lnrscale),guess_lnpscale/np.log(10.),\
         'r_scale,c=',r_scale,MYHALO.halo_concentration[i_snap],\
         'freeslope: c,alpha=',MYHALO.halo_concentration_freeslope[i_snap],\
         MYHALO.halo_innerslope_freeslope[i_snap])
        
        ## for more detailed analysis later, can just save the entire profile
        ##  (caution though -- the central profile can depend sensitively on centering, 
        ##    which ISN'T perfect here)
        rmin=r_s[is_dm_s][0]; rmax=r_s[is_dm_s][-1]; dlnr=np.log(rmax/rmin)/MYHALO.N_profile_radii;
        ln_rg=np.arange( np.log(rmin), np.log(rmax), dlnr )
        if(ln_rg.size<MYHALO.N_profile_radii): ln_rg=np.concatenate((ln_rg,np.array([ln_rg[-1]+dlnr])))
        if(ln_rg.size>MYHALO.N_profile_radii): ln_rg=ln_rg[0:np.int(MYHALO.N_profile_radii)]
        mg=np.interp( np.exp(ln_rg), r_s[is_dm_s], m_enc_dm )
        MYHALO.dm_mass_profile_rmin[i_snap] = rmin
        MYHALO.dm_mass_profile_rmax[i_snap] = rmax
        MYHALO.dm_mass_profile[i_snap,:] = mg
        
        ## gas-specific brick
        P = gadget.readsnap(snapdir,snum,0,\
            h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh)
        pos_g=P['p']; vel_g=P['v'];
        for j in [0,1,2]: pos_g[:,j] -= MYHALO.center[i_snap,j]
        for j in [0,1,2]: vel_g[:,j] -= MYHALO.com_vel[i_snap,j]
        r2_g = quadratic_sum(pos_g); 
        r_g = np.sqrt(r2_g);
        
        ## particles inside annuli, for inflow/outflow calculation:
        r_shells = r_200 * MYHALO.fraction_of_rvir_for_inflow_shells
        for i,r_shell in zip(range(r_shells.size),r_shells):
            shell_thick = 0.1
            ok_g = (r_g < r_shell/(1.+shell_thick)) & (r_g < r_shell*(1.+shell_thick))
            m_g=P['m'][ok_g]; h_g=P['h'][ok_g]; rho_g=P['rho'][ok_g]; zm_g=P['z'][ok_g,0];
            temp_g=gadget.gas_temperature(P['u'][ok_g],P['ne'][ok_g]);
            
            pos_g_ok=pos_g[ok_g,:]; vel_g_ok=vel_g[ok_g,:]; r_g_ok=r_g[ok_g]; v_r=0.*r_g_ok;
            for j in [0,1,2]: v_r += pos_g_ok[:,j]*vel_g_ok[:,j] / r_g_ok;
            h2=h_g*h_g; wt_h2_inv = 1./np.sum(h2)
            mdot_annulus_i = 4.*np.pi*r_shell*r_shell*wt_h2_inv * rho_g*h2*v_r
            delta_r = r_shell*(1.+shell_thick - 1./(1.+shell_thick))
            mdot_annulus_i = m_g * v_r / delta_r
                # delta_A = 4pi*h2*wt_h2_inv
                # rho __should be__ (since compressing to thin layer and don't want to double-count)
                ##  rho = dm/dV = dm/(dA*dr) ~ m_i / (delta_A * delta_r)
                # mdot = dm/(dr*dA) * dA*v_r ~ m_i * v_r/delta_r
            pdot_annulus_i = mdot_annulus_i * v_r
            mdot_metal_annulus_i = zm_g * mdot_annulus_i
            
            ## mass flux
            MYHALO.mdot_inflow[i_snap,i] = -np.sum(mdot_annulus_i[v_r < 0.])
            MYHALO.mdot_outflow[i_snap,i] = np.sum(mdot_annulus_i[v_r > 0.])
            ## momentum flux
            MYHALO.pdot_inflow[i_snap,i] = -np.sum(pdot_annulus_i[v_r < 0.])
            MYHALO.pdot_outflow[i_snap,i] = np.sum(pdot_annulus_i[v_r > 0.])
            ## metal flux
            MYHALO.mdot_inflow_metal[i_snap,i] = -np.sum(mdot_metal_annulus_i[v_r < 0.])
            MYHALO.mdot_outflow_metal[i_snap,i] = np.sum(mdot_metal_annulus_i[v_r > 0.])

            hot_thold_t = hot_gas_virial_threshold * temp_virial
            ## mass flux in 'hot' component 
            MYHALO.mdot_inflow_hot[i_snap,i] = -np.sum(mdot_annulus_i[(v_r < 0.) & (temp_g > hot_thold_t)])
            MYHALO.mdot_outflow_hot[i_snap,i] = np.sum(mdot_annulus_i[(v_r > 0.) & (temp_g > hot_thold_t)])
            ## momentum flux in 'hot' component 
            MYHALO.pdot_inflow_hot[i_snap,i] = -np.sum(pdot_annulus_i[(v_r < 0.) & (temp_g > hot_thold_t)])
            MYHALO.pdot_outflow_hot[i_snap,i] = np.sum(pdot_annulus_i[(v_r > 0.) & (temp_g > hot_thold_t)])
            ## metal flux in 'hot' component 
            MYHALO.mdot_inflow_metal_hot[i_snap,i] = -np.sum(mdot_metal_annulus_i[(v_r < 0.) & (temp_g > hot_thold_t)])
            MYHALO.mdot_outflow_metal_hot[i_snap,i] = np.sum(mdot_metal_annulus_i[(v_r > 0.) & (temp_g > hot_thold_t)])

            print('r in/out = ',r_shell)
            print(' .. mdot in/out ',MYHALO.mdot_inflow[i_snap,i],MYHALO.mdot_outflow[i_snap,i],\
                '.. metal ..',MYHALO.mdot_inflow_metal[i_snap,i],MYHALO.mdot_outflow_metal[i_snap,i])
            print(' .. mhot in/out ',MYHALO.mdot_inflow_hot[i_snap,i],MYHALO.mdot_outflow_hot[i_snap,i],\
                '.. metal ..',MYHALO.mdot_inflow_metal_hot[i_snap,i],MYHALO.mdot_outflow_metal_hot[i_snap,i])

        ## all gas inside halo:
        ok_g = (r_g < r_200)
        m_g=P['m'][ok_g]; id_g=P['id'][ok_g]; zm_g=P['z'][ok_g,0]; r_g=r_g[ok_g]
        temp_g=gadget.gas_temperature(P['u'][ok_g],P['ne'][ok_g]);

        ## total gas mass, and mass in various sub-categories
        MYHALO.m_gas_halo[i_snap] = np.sum(m_g)
        print('log10(m_gas) = ',np.log10(MYHALO.m_gas_halo[i_snap]))
        for i,temp_cut in zip(range(MYHALO.temperature_threshold_list.size),MYHALO.temperature_threshold_list):
            ok = (temp_g < temp_cut)
            MYHALO.m_gas_halo_tempcut[i_snap,i] = np.sum(m_g[ok])
            MYHALO.metal_gas_halo_tempcut[i_snap,i] = np.sum(m_g[ok]*zm_g[ok])
            print('log10(mgas<T)=',np.log10(MYHALO.m_gas_halo_tempcut[i_snap,i]),' T=',temp_cut)
            print('log10(Z<T)=',np.log10(MYHALO.metal_gas_halo_tempcut[i_snap,i]/0.02/MYHALO.m_gas_halo_tempcut[i_snap,i]),' T=',temp_cut)

        ## gas mass profile in halo
        s=np.argsort(r_g); rg_s=r_g[s]; mg_s=m_g[s]; mg_t=np.cumsum(mg_s);
        rmin=np.interp( 0.025*mg_t[-1], mg_t, rg_s );
        rmax=0.99*np.max(rg_s); r_eff=np.interp( 0.5*mg_t[-1], mg_t, rg_s );
        dlnr=0.02*np.log(rmax/rmin); ln_rg=np.arange( np.log(rmin), np.log(rmax), dlnr )
        ln_mg=np.interp( ln_rg, np.log(rg_s), np.log(mg_t) )
        sigma_g=np.exp(ln_mg)/(2.*np.pi*np.exp(2.*ln_rg)*2./3.) * \
            np.abs(my_simple_derivative( ln_mg, ln_rg ));
        ok=(sigma_g > 0.) & (np.isnan(sigma_g)==False);

        ##pylab.plot(ln_rg[ok]/np.log(10.),np.log10(sigma_g[ok]),'o--')
        
        ## for more detailed analysis later, can just save the entire profile
        ##  (caution though -- the central profile can depend sensitively on centering )
        rmin=rg_s[0]; rmax=rg_s[-1]; dlnr=np.log(rmax/rmin)/MYHALO.N_profile_radii;
        ln_rg=np.arange( np.log(rmin), np.log(rmax), dlnr )
        if(ln_rg.size<MYHALO.N_profile_radii): ln_rg=np.concatenate((ln_rg,np.array([ln_rg[-1]+dlnr])))
        if(ln_rg.size>MYHALO.N_profile_radii): ln_rg=ln_rg[0:np.int(MYHALO.N_profile_radii)]
        mg=np.interp( np.exp(ln_rg), rg_s, mg_t )
        MYHALO.gas_mass_profile_rmin[i_snap] = rmin
        MYHALO.gas_mass_profile_rmax[i_snap] = rmax
        MYHALO.gas_mass_profile[i_snap,:] = mg
            
            
        ## star-specific brick
        P = gadget.readsnap(snapdir,snum,4,\
            h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh)
        ## first check if there are any stars!
        if (P['k'] < 1):
            id_baryons=id_g
            m_baryons=m_g
            zm_baryons=zm_g
        else:
            if (P['m'].size<=0):
                id_baryons=id_g
                m_baryons=m_g
                zm_baryons=zm_g
            else:
                pos_s=P['p'];
                for j in [0,1,2]: pos_s[:,j] -= MYHALO.center[i_snap,j]
                r2_s = quadratic_sum(pos_s); 
                r_s = np.sqrt(r2_s);
        
                ## all stars inside halo:
                ok_s = (r_s < r_200)
                m_s=P['m'][ok_s]; id_s=P['id'][ok_s]; aform_s=P['age'][ok_s]; zm_s=P['z'][ok_s,0]
                pos_s=pos_s[ok_s,:]; vel_s=P['v'][ok_s,:]; r_s=r_s[ok_s];
                for j in [0,1,2]: vel_s[:,j] -= MYHALO.com_vel[i_snap,j]
            
                ## total stellar mass
                m_s_total = np.sum(m_s)
                MYHALO.m_stars_halo[i_snap] = m_s_total
                MYHALO.metal_stars_halo[i_snap] = np.sum(m_s*zm_s)
                ## get mass in new stars formed since last snapshot
                if (i_snap > 0):
                    a_previous_snapshot = MYHALO.a_scale[i_snap-1]
                else:
                    a_previous_snapshot = 0.
                ok = (aform_s > a_previous_snapshot)
                MYHALO.m_stars_new[i_snap] = np.sum(m_s[ok])
                MYHALO.metal_stars_new[i_snap] = np.sum(m_s[ok]*zm_s[ok])
                print('Total/new log10(m_stars) = ',np.log10(MYHALO.m_stars_halo[i_snap]),\
                    np.log10(MYHALO.m_stars_new[i_snap]))
                print('Total/new star metallicity: ',MYHALO.metal_stars_halo[i_snap],\
                    MYHALO.metal_stars_new[i_snap])
                
                ## make full baryon id list
                id_baryons = np.concatenate((id_g,id_s))
                m_baryons = np.concatenate((m_g,m_s))
                zm_baryons = np.concatenate((zm_g,zm_s))

                ## get angular momentum vectors for each star
                j_spec = np.cross( pos_s, vel_s , axis=1 )
                for j in [0,1,2]: MYHALO.j_stars[i_snap,j] = np.sum( m_s*j_spec[:,j] )
            
                ## get stellar effective radius and mass profile
                if (r_s.size > 10):
                    s=np.argsort(r_s); rs_s=r_s[s]; ms_s=m_s[s]; ms_t=np.cumsum(ms_s); fs_t=ms_t/ms_t[-1];
                    rmin=np.interp( 0.1*ms_t[-1], ms_t, rs_s );
                    rmax=np.interp( 0.9*ms_t[-1], ms_t, rs_s );
                    MYHALO.stellar_reff[i_snap] = np.interp( 0.5*ms_t[-1], ms_t, rs_s );
                    r_eff = MYHALO.stellar_reff[i_snap]
                    #
                    #dlnr=0.02*np.log(rmax/rmin); ln_rg=np.arange( np.log(rmin), np.log(rmax), dlnr )
                    #ln_mg=np.interp( ln_rg, np.log(rs_s), np.log(ms_t) )
                    #sigma_g=np.exp(ln_mg)/(2.*np.pi*np.exp(2.*ln_rg)*2./3.) * \
                    #    np.abs(my_simple_derivative( ln_mg, ln_rg ));
                    #p00 = [np.interp(r_eff,np.exp(ln_rg),np.log10(sigma_g)), np.log10(r_eff), 1.]
                    #MYHALO.stellar_sersic_parameters[i_snap,:],pcov = \
                    #    scipy.optimize.curve_fit(fitfun.sersic,ln_rg/np.log(10.),np.log10(sigma_g),p0=p00);
                    #q0 = p00;
                    #p00 = [q0[0]-1.,q0[1]+0.5,1. , q0[0]+1.,q0[1]-0.5,2.]
                    #MYHALO.stellar_sersicdisk_parameters[i_snap,:],pcov = \
                    #    scipy.optimize.curve_fit(fitfun.sersic_disk,ln_rg/np.log(10.),np.log10(sigma_g),p0=p00);
                    #
                    #print 'sersic fit ',MYHALO.stellar_sersic_parameters[i_snap,:]
                    #print 'sersic-disk fit ',MYHALO.stellar_sersicdisk_parameters[i_snap,:]
                    #pylab.plot(ln_rg/np.log(10.),np.log10(sigma_g),'o--')
                    #q=MYHALO.stellar_sersic_parameters[i_snap,:]
                    #pylab.plot(ln_rg/np.log(10.),fitfun.sersic(ln_rg/np.log(10.),q[0],q[1],q[2]),'r:')
                    #q=MYHALO.stellar_sersicdisk_parameters[i_snap,:]
                    #pylab.plot(ln_rg/np.log(10.),fitfun.sersic_disk(ln_rg/np.log(10.),q[0],q[1],q[2],q[3],q[4],q[5]),'m-')
                    # 
                    # -- sersic fitting can be done in post; it's rather unstable to do here --
                    ## for more detailed analysis later, can just save the entire profile
                    ##  (caution though -- the central profile can depend sensitively on centering )
                    rmin=rs_s[0]; rmax=rs_s[-1]; dlnr=np.log(rmax/rmin)/MYHALO.N_profile_radii;
                    ln_rg=np.arange( np.log(rmin), np.log(rmax), dlnr )
                    if(ln_rg.size<MYHALO.N_profile_radii): ln_rg=np.concatenate((ln_rg,np.array([ln_rg[-1]+dlnr])))
                    if(ln_rg.size>MYHALO.N_profile_radii): ln_rg=ln_rg[0:np.int(MYHALO.N_profile_radii)]
                    mg=np.interp( np.exp(ln_rg), rs_s, ms_t )
                    MYHALO.stellar_mass_profile_rmin[i_snap] = rmin
                    MYHALO.stellar_mass_profile_rmax[i_snap] = rmax
                    MYHALO.stellar_mass_profile[i_snap,:] = mg

                else:
                    MYHALO.stellar_reff[i_snap] = np.median(r_s)


        ## now some comparisons of our baryon list with baryons previously in halo
        sub_id_now, sub_id_prev, have_now_not_in_prev, have_prev_not_in_now = \
            mm.compile_matched_ids( id_baryons, id_baryons_in_previous_halo_timestep )
        # ok used it so can reset the 'previous timestep' vector to the current timestep
        id_baryons_in_previous_halo_timestep = id_baryons
        
        # total mass in material not in halo previous timestep 
        #  (from this, and gas/star mass, can reconstruct total mass that was in both, 
        #    and total mass lost from halo in previous timestep, though we also have that here)
        MYHALO.m_baryon_new_this_timestep[i_snap] = np.sum(m_baryons[have_now_not_in_prev])
        MYHALO.metals_baryon_new_this_timestep[i_snap] = np.sum(zm_baryons[have_now_not_in_prev] * \
            m_baryons[have_now_not_in_prev])

        # the material in the previous timestep is definitely already 'old news': 
        #   what about material which is new, but has been in the halo before?
        id_new_baryons = id_baryons[have_now_not_in_prev]
        m_new_baryons = m_baryons[have_now_not_in_prev]
        metals_new_baryons = zm_baryons[have_now_not_in_prev]
        
        sub_id_now, sub_id_prev, have_now_not_in_any, have_any_not_in_now = \
            mm.compile_matched_ids( id_new_baryons, id_baryons_in_any_halo_timestep )
        MYHALO.m_baryon_totally_new[i_snap] = np.sum(m_new_baryons[have_now_not_in_any])
        MYHALO.metals_baryon_totally_new[i_snap] = np.sum(metals_new_baryons[have_now_not_in_any] * \
            m_new_baryons[have_now_not_in_any])
        
        # add these (really) new particles to the list that have been in the halo before
        id_baryons_in_any_halo_timestep = np.concatenate((id_baryons_in_any_halo_timestep, \
            id_new_baryons[have_now_not_in_any]))
        
        print('mass/metal new to halo=',MYHALO.m_baryon_totally_new[i_snap],MYHALO.metals_baryon_totally_new[i_snap])
        print(' ')
        print(' ')
        sys.stdout.flush()


    #for key in MYHALO.__dict__.keys():
    #    print key
    #    print MYHALO.__dict__[key]

    ## open the file to record all saved properties
    fname = snapdir+'/halo_prop_vtime.hdf5'
    print('writing to file : ',fname)
    outfi = h5py.File(fname,'w')
    for key in MYHALO.__dict__.keys():
        dataset = outfi.create_dataset(key,data=MYHALO.__dict__[key])
    outfi.close()
    
    ## and make sure we dump the list of IDs that have appeared in the halo before, 
    ##   so that we can use it to restart the routine (if e.g. the snapshots are split 
    ##   across different machines or can't be loaded all at once)
    fname = snapdir+'/halo_prop_vtime_ids.hdf5'
    outfi_id = h5py.File(fname,'w')
    dataset = outfi_id.create_dataset('id_baryons_in_previous_halo_timestep',data=id_baryons_in_previous_halo_timestep)
    dataset = outfi_id.create_dataset('id_baryons_in_any_halo_timestep',data=id_baryons_in_any_halo_timestep)
    outfi_id.close()    
    print('done!')
    


def load_pos_vel(sdir,snum,ptypes=[0,1,2,3,4,5], 
        cosmological=1,skip_bh=1,four_char=0,use_h0=1):

	ma=np.zeros([0],dtype=float);
	ta=np.zeros([0],dtype=float);
	pa=np.zeros([0,3],dtype=float);
	va=np.zeros([0,3],dtype=float);
	for ptype in ptypes:
		P=gadget.readsnap(sdir,snum,ptype,\
		    h0=use_h0,four_char=four_char,cosmological=cosmological,skip_bh=skip_bh);
		if (P['k']==1): 
		    if(np.array(P['m']).size>0):
		        ma = np.concatenate((ma,P['m']))
		        ta = np.concatenate((ta,np.zeros((P['m'].size))+ptype))
		        pa = np.concatenate((pa,P['p']))
		        va = np.concatenate((va,P['v']))
	return ma,ta,pa,va


def quadratic_sum(x):
    return x[:,0]*x[:,0] + x[:,1]*x[:,1] + x[:,2]*x[:,2]


def nfw(log_r, log_norm, log_rscale):
        xc = log_r - log_rscale;
        return log_norm - xc - 2.*np.log(1.+np.exp(xc));


def nfw_freeinnerslope(log_r, log_norm, log_rscale, slope):
        xc = log_r - log_rscale;
        return log_norm - slope*xc - (3.-slope)*np.log(1.+np.exp(xc));


def tvir_of_halo( mhalo_code, rhalo_code ):
    msun = 1.989e33
    kpc = 3.086e21
    ggrav = 6.67e-8
    meanweight = 0.6
    mp = 1.67e-24
    kB = 1.3807e-16

    mhalo = mhalo_code*1.0e10 * msun
    rhalo = rhalo_code * kpc
    vc2 = ggrav*mhalo/rhalo 
    tvir = meanweight*mp * vc2 / (2.*kB)
    
    return tvir
    
    
def interp_w_extrap(x_new, x_old, y_old):
    y = np.interp(x_new,x_old,y_old)
    N = x_old.size - 1
    out = (x_new < x_old[0])
    if (y[out].size > 0):
        y[out] = y_old[0] + (y_old[1]-y_old[0])/(x_old[1]-x_old[0]) * (x_new[out]-x_old[0])
    out = (x_new > x_old[N])
    if (y[out].size > 0):
        y[out] = y_old[N] + (y_old[N]-y_old[N-1])/(x_old[N]-x_old[N-1]) * (x_new[out]-x_old[N])
    return y
    

def my_simple_derivative(y, x):
    xm=x[0:x.size-1]
    xp=x[1:x.size]
    x_mid = 0.5*(xm+xp)
    dy_dx = np.diff(y)/np.diff(x)
    dy_dx_xin = interp_w_extrap( x, x_mid, dy_dx )
    return dy_dx_xin



def look_at_halo():
    import matplotlib
    import matplotlib.pylab as pylab
    import matplotlib.pyplot as plot
    pylab.close('all')
    import warnings
    warnings.filterwarnings('ignore')
    
    sdir_master='/Users/phopkins/Documents/work/plots/zooms'

    SHOW_SFR_MAIN_SEQUENCE = 0
    SHOW_SSFR_REDSHIFT = 0
    SHOW_KSLAW = 1
    SHOW_KSLAW_NOFB = 0

    sdirs=['m12_mr']
    sdirs=['m11_hr']

    if ((SHOW_KSLAW)|(SHOW_KSLAW_NOFB)):
        fname_fig = 'zoom_kslaw.pdf'
        if (SHOW_KSLAW_NOFB): 
            fname_fig = 'zoom_kslaw_nofb.pdf'
        plot.figure(1,figsize=(8.,6.))
        charsi=20.
        matplotlib.rcParams.update({'font.size':18})
        sdirs=['hires_bh', 'm10_dw', 'm12_mr', 'm11_hr', 'm11_lr','m13_mr','m12_lr' ]

        #sdirs=[default_B1,  default_m10, default_m12,      default_m11,      default_m09,      default_m13,  default_m12v]#,default_m14]
        sdirs=['B1_Sep10_mhd', 'm10_lr_Jan2014', 'm12qq_Jan2014', 'm11_lr_Jan2014', 'm14_Jan2014',  'm13_mr','m12v_3_Jan2014']#,'m14_Jan2014' ]
        sdirs=['B1_11_Jan2014', 'm10_lr_Jan2014', 'm12qq_Jan2014','m11_hr_Jan2014','m13_mr','m13_mr_Jan2014','m12v_3_Jan2014','m14_Jan2014' ]
        sdirs=['B1_11_Jan2014', 'm10_hr_Jan2014', 'm12qq_Jan2014','m11_hr_Jan2014','m13_mr','m14_Jan2014','m12v_3_Jan2014']#,'m14_Jan2014' ]

        sdirs=['B1_11_Jan2014', 'm10_hr_Jan2014', 'm12v_3_Jan2014','m11_hr_Jan2014','m13_mr','m14_Jan2014','m12qq_Jan2014']#,'m14_Jan2014' ]
        #sdirs=['m13_mr_Jan2014','m13_mr','m13_mr_Aug8']

        if (SHOW_KSLAW_NOFB):
            sdirs=['m10_norad', 'm11_nofb', 'm13_nofb', 'z10n192_B1hr_cwef1en25', 'z10n192_B1hr_cwef2en50']
            sdirs=['hires_bh_nosn','m10_norad','hires_bh_nosn','m10_norad']
        marker_syms=['o','s','p','D','D','*','x','x','D']
        colors_vec=['k','b','r','g','g','c','m','m','y']

        marker_syms=['o','s','p','^','*','x','D','+','h','1','v','2','<','3','>','4','H']
        colors_vec=['black','blue','red','green','deepskyblue','darkviolet','orange','magenta','gold','sienna','pink','forestgreen','darkgray']

        marker_syms=['o','s','D','^','*','x','p','+','h','1','v','2','<','3','>','4','H']
        colors_vec=['black','blue','orange','green','deepskyblue','darkviolet','red','magenta','gold','sienna','pink','forestgreen','darkgray']


        x=np.arange(8.,10.5,0.1)
        x0=np.array([9.5,10.,10.5])
        lthick=3.
        pylab.axis([-0.6,4.,-5.8,2.])
        #pylab.axis([-0.6,5.,-5.8,4.])
        pylab.xlabel(r'$\log(\ \Sigma_{\rm gas}\ )\ \ [{\rm M}_{\odot}\,{\rm pc}^{-2}]$')
        pylab.ylabel(r'$\log(\ \Sigma_{\rm SFR}\ )\ \ [{\rm M}_{\odot}\,{\rm yr}^{-1}\,{\rm kpc}^{-2}]$')
        xx=np.arange(-4.,6.,0.1)
        #pylab.plot(xx,-1.+(xx-2.)*1.5,'k--',linewidth=lthick)
        #pylab.plot(xx,-1.+(xx-2.)*1.75,'k-',linewidth=lthick)
        #pylab.plot(xx,-1.+(xx-2.)*2.0,'k--',linewidth=lthick)

        xx=10.**np.arange(5.,13.,0.1)
        y1=np.log10((1.e8/0.41e9)*((xx/1.e8)**(1.71)))
        y1=np.log10((1.e8/0.21e9)*((xx/1.e8)**(1.76)))
        y2=np.log10((1.e8/0.63e9)*((xx/1.e8)**(1.40)))
        y3=np.log10((1.e8/3.00e9)*((xx/1.e8)**(1.91)))
        yp=y1; yp[y2>yp]=y2[y2>yp]; yp[y3>yp]=y3[y3>yp]; 
        ym=y1; ym[y2<ym]=y2[y2<ym]; ym[y3<ym]=y3[y3<ym]; 
        x=np.log10(xx)-6.; y0=0.5*(yp+ym); dy=0.5*(yp-ym);
        dyyp=np.sqrt(dy*dy+0.6**2.); dyym=np.sqrt(dy*dy+0.6**2.); yp=y0+dyyp; ym=y0-dyym;
        #pylab.fill_between(x,ym,yp,facecolor='pink',alpha=0.5)

        log_surfacedensity=[\
        0.506100, 0.745300, 0.905100,  1.08900,  1.19400,  1.45700,  1.77000,  2.10700,  2.41300,  2.58500,  2.84300,\
         3.17400,  3.71300,  4.18500,  3.95900,  3.57900,  3.19300,  2.70900,  2.26800,  1.91800,  1.33000,  1.12800,\
        0.656000, 0.404400, 0.183300,  0.09706, -0.24000]
        log_surfacesfr=[\
        -4.49700, -4.11300, -3.54200, -3.05700, -2.53300, -2.22500, -1.88000, -1.54200, -1.09600, -0.64930, -0.11810, 0.350900, 0.973400,  1.39500,\
         2.22200,  1.62200,  1.02200, 0.800600, 0.231700, -0.62310, -1.04400, -1.39100, -2.03700, -2.63000, -3.24600, -3.75500, -4.17000]

        log_surfacedensity=[\
        0.506100, 0.745300, 0.905100,  1.08900,  1.19400,  1.45700,  1.77000,  2.10700,  2.41300,  2.58500,  2.84300,\
         3.17400,  3.71300,  4.18500,  3.95900,  3.57900,  3.19300,  2.70900,  2.26800,  1.91800,  1.33000,  1.12800,\
        0.656000, 0.404400, 0.183300,  0.09706, -0.24000]
        log_surfacesfr=[\
        -4.49700, -4.11300, -3.54200, -3.05700, -2.53300, -2.22500, -1.88000, -1.54200, -1.09600, -0.64930,  0.36310,\
        0.850900, 1.373400,  1.89500,  2.82200,  2.32200,  1.62200, 1.200600, 0.431700, -0.32310, -1.04400, -1.39100,\
        -2.03700, -2.63000, -3.24600, -3.75500, -4.17000]
        #pylab.plot(np.array(log_surfacedensity)+0.3,np.array(log_surfacesfr),'k-',linewidth=2.)

        x=np.arange(-1.,4.5,0.1)
        log_surfacedensity=np.array([\
        -1., 0.0, 0.506100, 0.745300, 0.905100,  1.08900,  1.19400,  1.45700,  1.77000,  2.10700,  2.41300,  2.58500,  2.84300,\
         3.17400,  3.71300])
        log_surfacesfr=np.array([\
        -6., -6.0, -4.49700, -4.11300, -3.54200, -3.05700, -2.53300, -2.22500, -1.88000, -1.54200, -1.09600, -0.64930,  0.36310,\
        0.850900, 1.373400])
        s=np.argsort(log_surfacedensity); log_surfacedensity=log_surfacedensity[s]; log_surfacesfr=log_surfacesfr[s];
        y1=np.interp(x,log_surfacedensity,log_surfacesfr)
        #pylab.plot(np.array(log_surfacedensity)+0.3,np.array(log_surfacesfr),'k-',linewidth=2.)
        log_surfacedensity=np.array([\
         4.18500,  3.95900,  3.57900,  3.19300,  2.70900,  2.26800,  1.91800,  1.33000,  1.12800,\
        0.656000, 0.404400, 0.183300,  0.09706, -0.24000, -0.75, -1.])
        log_surfacesfr=np.array([\
        1.89500,  2.82200,  2.32200,  1.62200, 1.200600, 0.431700, -0.32310, -1.04400, -1.39100,\
        -2.03700, -2.63000, -3.24600, -3.75500, -4.17000, -5.30, -6.])
        s=np.argsort(log_surfacedensity); log_surfacedensity=log_surfacedensity[s]; log_surfacesfr=log_surfacesfr[s];
        y2=np.interp(x,log_surfacedensity,log_surfacesfr)
        #pylab.plot(np.array(log_surfacedensity)+0.3,np.array(log_surfacesfr),'k-',linewidth=2.)
        pylab.fill_between(x+0.3,y1,y2,facecolor='yellow',alpha=0.5)


        log_surfacedensity=[\
          0.392300,  0.531600,  0.726100,  0.755900,  0.779300,  0.845900,\
          0.955500,   1.15700,   1.22300,   1.31400,   1.52200,   1.69800,\
           1.84400,   2.12500,   2.32600,   2.47300,   2.65000,   2.73500,\
           2.73900,   2.81700,   2.98000,   3.53700,   4.12300]
        log_surfacesfr=[\
          -5.01800,  -4.63200,  -4.23200,  -4.03100,  -3.79900,  -3.60700,\
          -3.41400,  -3.25300,  -2.96800,  -2.69000,  -2.53600,  -2.28200,\
          -1.98200,  -1.60500,  -1.46700,  -1.22800,  -1.09000, -0.897800,\
         -0.473100, -0.164600,  0.390600,  0.812700,   1.35000]
        #pylab.plot(np.array(log_surfacedensity)+0.2,np.array(log_surfacesfr),'r-',linewidth=4.)

        log_surfacedensity=[\
            1.75700, 1.86100, 1.87200, 2.43900, 2.56300, 2.56500, 2.61500, 2.65100, 2.67000, 2.69600, 2.76900,\
            3.02100, 3.02100, 3.21400, 3.22800, 3.23200, 3.31300, 3.36000, 3.38700, 3.44500, 3.49300, 3.54700,\
            3.55800, 3.62200, 3.63900, 3.75500, 3.75900, 3.89100, 3.91300, 4.01300, 4.30600]
        log_surfacesfr=[\
          -0.099240, 0.24600, -0.5681, 0.30570, 0.26690, 0.63550, 0.75060, 0.71980, 0.950100, 1.13400, 1.11100,\
            1.21700, 1.17900, 1.40100, 1.44000, 1.15500, 1.54700, 1.29300, 1.47000, 1.47000, 2.16800, 1.94500,\
            1.68400, 1.75300, 1.55300, 1.54500, 2.32100, 2.19000, 2.19800, 2.29700, 2.54200]
        #pylab.plot(np.array(log_surfacedensity)+0.3,np.array(log_surfacesfr),'o',color='red',linewidth=4.)


    if (SHOW_SFR_MAIN_SEQUENCE):
        fname_fig = 'zoom_sfr_mainsequence.pdf'
        plot.figure(1,figsize=(8.,6.))
        charsi=20.
        matplotlib.rcParams.update({'font.size':18})
        sdirs=['hires_bh', 'm10_dw', 'm12_mr', 'm11_hr', 'm11_lr','m13_mr','m12_lr' ]
        #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_hr', 'm09_dw', 'm13_gas_lr','m13_gas_lr','m14_tst']
        marker_syms=['o','s','p','D','D','*','x','x','D']
        colors_vec=['k','b','r','g','g','c','m','m','y']
        
        x=np.arange(8.,10.5,0.1)
        x0=np.array([9.5,10.,10.5])
        lthick=3.
        yscatter=np.log10(2.0)

        ## z = 0
        ssfr=np.array([7.e-11,5.e-11,3.e-11]); y=np.interp(x,x0,np.log10(ssfr)) + x
        pylab.plot(x,y,'k-',linewidth=lthick)
        #pylab.plot(x,y-yscatter,'k--',linewidth=lthick)
        #pylab.plot(x,y+yscatter,'k--',linewidth=lthick)
        ## z = 1
        ssfr=np.array([9.6e-10,8.e-10,4.7e-10]); y=np.interp(x,x0,np.log10(ssfr))+x
        pylab.plot(x,y,'b-',linewidth=lthick)
        #pylab.plot(x,y-yscatter,'b--',linewidth=lthick)
        #pylab.plot(x,y+yscatter,'b--',linewidth=lthick)
        ## z = 2
        ssfr=np.array([2.1e-9,1.8e-9,1.4e-9]); y=np.interp(x,x0,np.log10(ssfr)) + x
        pylab.plot(x,y,'g-',linewidth=lthick)
        #pylab.plot(x,y-yscatter,'g--',linewidth=lthick)
        #pylab.plot(x,y+yscatter,'g--',linewidth=lthick)
        ## z = 4
        ssfr=np.array([4.5e-9,4.4e-9,4.0e-9]); y=np.interp(x,x0,np.log10(ssfr)) + x
        pylab.plot(x,y,'r-',linewidth=lthick)
        #pylab.plot(x,y-yscatter,'r--',linewidth=lthick)
        #pylab.plot(x,y+yscatter,'r--',linewidth=lthick)

    for i_sdir,sdir0 in zip(range(np.array(sdirs).size),sdirs):
    
        snapdir = sdir_master+'/'+sdir0
        N_snaps = 441 # need to know this ahead of time, will look for every snap from 0-this number
        hot_gas_virial_threshold = 0.5 # fraction of t_vir to qualify as 'hot' gas (for inflow/outflow)
        virial_threshold = 200. # standard virial overdensity threshold

        ## initialize variables needed below
        MH = halo(N_snaps)
        fname = snapdir+'/halo_prop_vtime.hdf5'
        print('getting data from file : ',fname)
        infi = h5py.File(fname,'r')
        #print infi.keys()
        for key in infi.keys():
            MH.__dict__[key] = np.copy(infi[key])
        infi.close()
        tmax = 14.

        msun = 1.989e33
        munit = 1.e10*msun
        Ggrav = 6.67e-8
        kpc = 3.086e21
        runit = kpc
        kms = 1.0e5
        j_unit = kms*runit
        yr = 3.155e7
        mdot_unit = munit/runit*kms / (msun/yr)

        z=1./MH.a_scale-1.
        t=MH.time
        #t=MH.a_scale
        v_200 = np.sqrt(Ggrav*MH.m_200*munit/(MH.r_200*runit)) / kms
        log_mh = np.log10(MH.m_200) + np.log10(munit/msun)
        log_mg = np.log10(MH.m_gas_halo) + np.log10(munit/msun)
        log_ms = np.log10(MH.m_stars_halo) + np.log10(munit/msun)
        log_ms_new = np.log10(MH.m_stars_new) + np.log10(munit/msun)
        sfr = 10.**(log_ms_new[1:log_ms_new.size]) / (np.diff(t)*1.e9)
        sfr_t = MH.time[0:MH.time.size-1]+0.5*np.diff(MH.time)
        ok=(np.isnan(sfr)==False) & (sfr>0)
        sfr_n = np.interp(MH.time,sfr_t[ok],sfr[ok]) / 0.7
        r_eff = 1.0*MH.stellar_reff
        if(SHOW_KSLAW):
            r_eff *= 1.5;
        if(SHOW_KSLAW_NOFB):
            if(i_sdir >= 2):
                r_eff *= 0.5;
        sigma_sfr = sfr_n / (np.pi*r_eff*r_eff)
        dlnr=np.log(MH.gas_mass_profile_rmax/MH.gas_mass_profile_rmin)/np.float(MH.N_profile_radii)
        mgas_enc_r_eff = 0.*r_eff
        for i in range(MH.gas_mass_profile.shape[0]):
            if ((MH.gas_mass_profile_rmin[i]>0) & (MH.gas_mass_profile_rmax[i]>0)):
                r=np.exp(np.arange(np.log(MH.gas_mass_profile_rmin[i]),np.log(MH.gas_mass_profile_rmax[i]),dlnr[i]))
                r*=np.sqrt(2./3.); r=r[0:np.int(MH.N_profile_radii)];
                ok=(r <= r_eff[i])
                if (MH.gas_mass_profile[i,ok].size > 0):
                    mgas_enc_r_eff[i] = np.max(MH.gas_mass_profile[i,ok])
                #r=r[0:np.int(MH.N_profile_radii)]
                #wt=MH.gas_mass_profile[i,:]*4.*np.pi*r*r*r*dlnr[i]
                #print MH.gas_mass_profile[i,:]
                #mgas_enc_r_eff[i] = np.sum(wt[(r <= r_eff[i])])
        sigma_gas_surface = 2. * mgas_enc_r_eff / (np.pi*r_eff*r_eff)
                
        #sigma_gas_surface = 
        Z_stars = np.log10(MH.metal_stars_halo/MH.m_stars_halo)
        Z_stars_new = np.log10(MH.metal_stars_new/MH.m_stars_new)
        log_m_new_baryons = np.log10(MH.m_baryon_new_this_timestep) + np.log10(munit/msun)
        Z_new_baryons = np.log10(MH.metals_baryon_new_this_timestep/MH.m_baryon_new_this_timestep)
        f_totally_new_baryons = MH.m_baryon_totally_new/MH.m_baryon_new_this_timestep
        Z_totally_new_baryons = np.log10(MH.metals_baryon_totally_new/MH.m_baryon_totally_new)
        vmax = np.sqrt(Ggrav*MH.v2_max_all*munit/runit) / kms
        vmax_dm = np.sqrt(Ggrav*MH.v2_max_dm*munit/runit) / kms
        jmag_total = np.sqrt(quadratic_sum(MH.j_total)) / (MH.m_200)
        jmag_dm = np.sqrt(quadratic_sum(MH.j_total_dm)) / (MH.m_200-MH.m_gas_halo-MH.m_stars_halo)
        j_dm = MH.j_total_dm
        j_baryons = MH.j_total - j_dm
        jmag_baryons = np.sqrt(quadratic_sum(j_baryons)) / (MH.m_gas_halo+MH.m_stars_halo)
        j_stars = MH.j_stars
        jmag_stars = np.sqrt(quadratic_sum(j_stars)) / (MH.m_stars_halo)
        spin = jmag_total*j_unit / (np.sqrt(2.)*MH.r_200*runit * v_200*kms)
        inflow_inner_outer = MH.mdot_inflow[:,0]/MH.mdot_inflow[:,3]
        outflow_inner_outer = MH.mdot_outflow[:,0]/MH.mdot_outflow[:,3]
        Z_outflow = np.log10(MH.mdot_outflow_metal/MH.mdot_outflow)
        Z_inflow = np.log10(MH.mdot_inflow_metal/MH.mdot_inflow)
        dmhalo_dt=0.*t; ok=(t>0.); dmhalo_dt[ok]=my_simple_derivative(10.**log_mh[ok],t[ok]*1.0e9)
        dmbaryon_dt=0.*t; ok=(t>0.); dmbaryon_dt[ok]=my_simple_derivative(10.**log_mg[ok],t[ok]*1.0e9)
        mdot_expected = dmhalo_dt * 0.162
    
        #pylab.plot(t,log_mh,'o-')
        #pylab.plot(t,np.log10(MH.r_200),'r-')
        #pylab.plot(log_mh,np.log10(MH.r_200),'ro-')
    
        """
        pylab.plot(t,vmax,'o')
        pylab.plot(t,vmax_dm,'.')
        """
        """
        pylab.axis([0.,tmax,0.,30.])
        pylab.plot(t,MH.halo_concentration,'o')
        pylab.plot(t,MH.halo_concentration_freeslope,'.')
        pylab.plot(t,MH.halo_concentration/MH.a_scale,'ro')
        """
        """
        pylab.axis([0.,tmax,-2.5,2.5])
        pylab.plot(t,MH.halo_innerslope_freeslope,'o')
        """
        """
        #t=MH.a_scale
        #pylab.axis([0.63,0.65,8.,10.])
        pylab.axis([8.,9.2,0.,10.])
        ok=(t > 0.63) & (t < 0.65)
        ii=np.arange(0,t.size,1)
        print ii[ok]
        pylab.plot(t,log_mh+np.log10(0.162),'ko')
        pylab.plot(t,log_mg,'b.')
        pylab.plot(t,log_ms,'g.')
        #pylab.plot(t,log_ms_new,'r.')
        ok=(np.isnan(sfr)==False) & (sfr>0)
        sfr_n = np.interp(t,sfr_t[ok],sfr[ok]) / 0.7
        pylab.plot(sfr_t[ok],sfr[ok]*500.,'r.')
        """
        
        if (SHOW_SFR_MAIN_SEQUENCE):
            ## need to do something slightly fancier to account for mass-loss here!
            ok=(z < 0.5) & (t>0)
            pylab.plot(log_ms[ok],np.log10(sfr_n[ok] * 0.7),'k.')
            ok=(z > 0.5) & (z < 1.5) & (t>0)
            ok=(z > 0.8) & (z < 1.5) & (t>0)
            pylab.plot(log_ms[ok],np.log10(sfr_n[ok]),'b.')
            ok=(z > 1.5) & (z < 2.5) & (t>0)
            ok=(z > 1.8) & (z < 2.5) & (t>0)
            pylab.plot(log_ms[ok],np.log10(sfr_n[ok]),'g.')
            ok=(z > 2.5) & (t>0)
            ok=(z > 3.5) & (z < 5.) & (t>0)
            pylab.plot(log_ms[ok],np.log10(sfr_n[ok]),'r.')

        if (SHOW_SSFR_REDSHIFT):
            ## need to do something slightly fancier to account for mass-loss here!
            ok=(z < 0.5) & (t>0)
            pylab.plot(log_ms[ok],np.log10(sfr_n[ok] * 0.7),'k.')
            ok=(z > 0.5) & (z < 1.5) & (t>0)
            ok=(z > 0.8) & (z < 1.5) & (t>0)
            pylab.plot(log_ms[ok],np.log10(sfr_n[ok]),'b.')
            ok=(z > 1.5) & (z < 2.5) & (t>0)
            ok=(z > 1.8) & (z < 2.5) & (t>0)
            pylab.plot(log_ms[ok],np.log10(sfr_n[ok]),'g.')
            ok=(z > 2.5) & (t>0)
            ok=(z > 3.5) & (z < 5.) & (t>0)
            pylab.plot(log_ms[ok],np.log10(sfr_n[ok]),'r.')
    
    
        if (SHOW_KSLAW)|(SHOW_KSLAW_NOFB):
            xx = np.log10(sigma_gas_surface)+10.-6.
            #if (sdir0=='m13_mr'): xx -= 2.
            #if (sdir0=='m10_dw'): xx -= 6.
            xoffset=1.7
            if(sdir0=='m14_Jan2014'): xoffset=0.6
            if(sdir0=='m10_hr_Jan2014'): xoffset=1.0
            if(sdir0=='m12qq_Jan2014'): xoffset=1.5
            if(sdir0=='B1_11_Jan2014'): xoffset=1.5
            if(sdir0=='m11_hr_Jan2014'): xoffset=1.2
            if(sdir0=='m12v_3_Jan2014'): xoffset=1.5
            if(sdir0=='m13_mr'): xoffset=1.5
            #sdirs=['m13_mr','m12v_3_Jan2014']
            ok=(np.random.rand(xx.size)<1.5)&(xx>=np.median(xx)-xoffset)
            pylab.scatter(xx[ok],np.log10(sigma_sfr[ok]),\
                marker=marker_syms[i_sdir],color=colors_vec[i_sdir],\
                facecolors='none')
            #print sigma_gas_surface
            #print mgas_enc_r_eff
            #print sigma_sfr
        
        """
        pylab.plot(t,Z_stars)
        pylab.plot(t,Z_stars_new,'r')
        """
        #pylab.plot(t,MH.stellar_reff)
        """
        pylab.plot(t,MH.center[:,0]/MH.a_scale,'k')
        pylab.plot(t,MH.center[:,1]/MH.a_scale,'b')
        pylab.plot(t,MH.center[:,2]/MH.a_scale,'r')
        """
        #pylab.plot(t,log_m_new_baryons)
        #pylab.plot(t,log_ms_new,'r')
    
        #pylab.plot(t,Z_new_baryons)
        #pylab.plot(t,Z_totally_new_baryons,'r')
    
        #pylab.plot(t,f_totally_new_baryons)
        """
        pylab.yscale('log')
        pylab.plot(t,jmag_total,'k')
        pylab.plot(t,jmag_dm,'b')
        pylab.plot(t,jmag_baryons,'g')
        pylab.plot(t,jmag_stars,'r')
        """
        #pylab.plot(t,spin)

        """
        pylab.yscale('log')
        #pylab.xscale('log') # if x-axis is log-a_scale, 'periodicity' looks more regular
        ok=(t>0)
        #pylab.axis([10.,11.,0.001,0.1])
        #pylab.axis([0.,9.,0.3,30.])
        pylab.plot(t[ok],MH.mdot_inflow[ok,0],'ko-')
        #pylab.plot(t[ok],MH.mdot_inflow_hot[ok,0],'ko--')
        pylab.plot(t[ok],MH.mdot_outflow[ok,0],'r.-')
        #pylab.plot(t[ok],MH.mdot_outflow_hot[ok,0],'r.--')
        pylab.plot(t[ok],MH.mdot_inflow[ok,3],'bo-')
        pylab.plot(t[ok],MH.mdot_outflow[ok,3],'g.-')
        """
        """
        pylab.yscale('log')
        ok=(t>0)
        #pylab.axis([0.5,3.,1.e-3,10.])
        pylab.plot(t[ok],MH.mdot_inflow[ok,0]*mdot_unit,'k,-')
        pylab.plot(t[ok],MH.mdot_outflow[ok,0]*mdot_unit,'r,-')
        #pylab.plot(t[ok],MH.mdot_outflow[ok,1]*mdot_unit,'m,-')
        #pylab.plot(t[ok],MH.mdot_inflow[ok,3]*mdot_unit*f_totally_new_baryons[ok],'mo-')
        #pylab.plot(t[ok],MH.mdot_inflow[ok,1]*mdot_unit*f_totally_new_baryons[ok],'mo-')
        #pylab.plot(t[ok],MH.mdot_inflow[ok,0]*mdot_unit,'mo-')
        pylab.plot(t[ok],mdot_expected[ok],'b.-')
        pylab.plot(t[ok],dmbaryon_dt[ok],'go-')
        #pylab.plot(t[ok],dmbaryon_dt[ok]*f_totally_new_baryons[ok],'ro-')
        """
        """
        pylab.yscale('log')
        ok=(t>0)
        pylab.axis([0.,tmax,0.01,10.])
        pylab.plot(t[ok],inflow_inner_outer[ok],'bo-')
        pylab.plot(t[ok],outflow_inner_outer[ok],'r.-')
        """
        """
        pylab.yscale('log')
        ok=(t>0)
        pylab.axis([0.,tmax,0.01,100.])
        pylab.plot(t[ok],0.*t[ok]+1.,'b--')
        for j,col in zip([0,3],['k-','r-','g--']):
            pylab.plot(t[ok],MH.mdot_inflow[ok,j]/MH.mdot_outflow[ok,j],col)
        """
        """
        pylab.yscale('log')
        ok=(t>0)
        pylab.plot(t[ok],-MH.pdot_inflow[ok,0],'ko-')
        #pylab.plot(t[ok],-MH.pdot_inflow_hot[ok,0],'ko--')
        pylab.plot(t[ok],MH.pdot_outflow[ok,0],'r.-')
        #pylab.plot(t[ok],MH.pdot_outflow_hot[ok,0],'r.--')
        pylab.plot(t[ok],-MH.pdot_inflow[ok,3],'bo-')
        pylab.plot(t[ok],MH.pdot_outflow[ok,3],'g.-')
        #pylab.plot(t[ok],10.**log_ms_new[ok],'m-')
        """
        """
        ok=(t>0)
        pylab.plot(t[ok],Z_inflow[ok,0],'ko-')
        pylab.plot(t[ok],Z_outflow[ok,0],'r.-')
        pylab.plot(t[ok],Z_inflow[ok,3],'ko--')
        pylab.plot(t[ok],Z_outflow[ok,3],'r.--')
        pylab.plot(t,Z_stars,'g--')
        pylab.plot(t,Z_stars_new,'b-')
        """
        """
        prof='stellar_mass'
        prof='gas_mass'
        #prof='dm_mass'
        #pylab.axis([0.,tmax,0.01,100.]); 
        m=MH.__dict__[prof+'_profile']; rmin=MH.__dict__[prof+'_profile_rmin']; rmax=MH.__dict__[prof+'_profile_rmax'];
        N=np.float(MH.N_profile_radii); dlnr=np.log(rmax/rmin)/N; 
        for j in range(rmin.size):
            if ((rmin[j]>0.) & (t[j]>13.5)):
                ln_rg=np.arange( np.log(rmin[j]), np.log(rmax[j]), dlnr[j] )
                if(ln_rg.size < N): ln_rg=np.concatenate((ln_rg,np.array([ln_rg[-1]+dlnr])))
                if(ln_rg.size > N): ln_rg=ln_rg[0:np.int(N)]
                ln_mg=np.log(m[j,:])
                sigma = np.exp(ln_mg)/(2.*np.pi*np.exp(2.*ln_rg)*2./3.) * \
                    np.abs(my_simple_derivative( ln_mg, ln_rg ));
                rho = np.exp(ln_mg)/(4.*np.pi*np.exp(3.*ln_rg)) * \
                    np.abs(my_simple_derivative( ln_mg, ln_rg ));

                pylab.plot(ln_rg/np.log(10.),np.log10(rho),'-')
        """
    
    pylab.savefig(fname_fig,transparent=True,bbox_inches='tight',pad_inches=0)




