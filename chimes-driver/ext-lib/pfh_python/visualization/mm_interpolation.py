import numpy as np
import gc
from mm_utilities import *


def interpolation_for_movie(dt, delta_t, \
	    m_i,h_i,c_i,zm_i,\
	    m_f,h_f,c_f,zm_f, \
	    x_in,y_in,z_in,vx_in,vy_in,vz_in, \
	    x_fn,y_fn,z_fn,vx_fn,vy_fn,vz_fn, \
	    s_i,s_f,nm_i_f,nm_f_i, \
	    cen_i,cen_f , use_polar_interpolation=1, 
	    polar_vs_cylindrical = True, 
	    set_dx_frame = 0., stars = 0):
    
    x_i = np.copy(x_in) - cen_i[0]
    y_i = np.copy(y_in) - cen_i[1]
    z_i = np.copy(z_in) - cen_i[2]
    vx_i = np.copy(vx_in) - (cen_f[0]-cen_i[0])/delta_t
    vy_i = np.copy(vy_in) - (cen_f[1]-cen_i[1])/delta_t
    vz_i = np.copy(vz_in) - (cen_f[2]-cen_i[2])/delta_t
    x_f = np.copy(x_fn) - cen_f[0]
    y_f = np.copy(y_fn) - cen_f[1]
    z_f = np.copy(z_fn) - cen_f[2]
    vx_f = np.copy(vx_fn) - (cen_f[0]-cen_i[0])/delta_t
    vy_f = np.copy(vy_fn) - (cen_f[1]-cen_i[1])/delta_t
    vz_f = np.copy(vz_fn) - (cen_f[2]-cen_i[2])/delta_t
    
    #for ui,uin,j in zip([x_i,y_i,z_i],[x_in,y_in,z_in],[0,1,2]): ui = uin - cen_i[j]
    #for uf,ufn,j in zip([x_f,y_f,z_f],[x_fn,y_fn,z_fn],[0,1,2]): uf = ufn - cen_f[j]
    #for vui,vuin,j in zip([vx_i,vy_i,vz_i],[vx_in,vy_in,vz_in],[0,1,2]): vui = vuin - (cen_f[j]-cen_i[j])/delta_t
    #for vuf,vufn,j in zip([vx_f,vy_f,vz_f],[vx_fn,vy_fn,vz_fn],[0,1,2]): vuf = vufn - (cen_f[j]-cen_i[j])/delta_t
    pos_interp_flag = np.array([0.])
    if(stars!=0): 
        r_i = np.sqrt(x_i[s_i]*x_i[s_i] + y_i[s_i]*y_i[s_i] + z_i[s_i]*z_i[s_i])
        r_f = np.sqrt(x_f[s_f]*x_f[s_f] + y_f[s_f]*y_f[s_f] + z_f[s_f]*z_f[s_f])
        d_r = np.sqrt((x_i[s_i]-x_f[s_f])**2 + (y_i[s_i]-y_f[s_f])**2 + (z_i[s_i]-z_f[s_f])**2)
        r_eff = 1. / (1./r_i + 1./r_f)
        dvx=(x_f[s_f]-x_i[s_i])/delta_t; dvy=(y_f[s_f]-y_i[s_i])/delta_t; dvz=(z_f[s_f]-z_i[s_i])/delta_t
        dv2_i=(vx_i[s_i]-dvx)**2 + (vy_i[s_i]-dvy)**2 + (vz_i[s_i]-dvz)**2
        dv2_f=(vx_f[s_f]-dvx)**2 + (vy_f[s_f]-dvy)**2 + (vz_f[s_f]-dvz)**2
        dv = np.sqrt(dv2_i + dv2_f); v0 = np.sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        pos_interp_flag = (dv/v0 > 1.0) & (dv > 100.) & (d_r/r_eff > 2.2) & (d_r > 0.01*set_dx_frame) & (r_eff < 0.2*set_dx_frame)

    if (use_polar_interpolation==1):
        # convert to cylindrical coordinates, with the axis defined as cumulative ang-mom axis inside each orbit
        r_i,vr_i,zc_i,vzc_i,phi_i,vphi_i,j_hat_i = xyz_to_cylindrical(x_i,y_i,z_i,vx_i,vy_i,vz_i)
        r_f,vr_f,zc_f,vzc_f,phi_f,vphi_f,j_hat_f = xyz_to_cylindrical(x_f,y_f,z_f,vx_f,vy_f,vz_f)
        if(stars!=0): 
            z0 = np.sqrt(zc_i[s_i]*zc_i[s_i] + zc_f[s_f]*zc_f[s_f])
            dz = np.abs(zc_i[s_i]-zc_f[s_f])
            #pos_interp_flag = pos_interp_flag & (z0/r_eff > 0.5) & (z0 > 0.5)
            #pos_interp_flag = pos_interp_flag & (dz/z0 > 1.0) & (r_eff > 0.5) & (r_eff < 3.) 
            pos_interp_flag = pos_interp_flag & (r_eff < 3.) & (dz > 0.1) & (dz/z0 > 1.) & (z0 > 0.1)
            
            
            #pos_interp_flag = pos_interp_flag & (z0/r_eff > 0.5) & (z0 > 0.25) & (dz/z0 > 0.75) & (r_eff < 3.)
        if(polar_vs_cylindrical == True):
            # convert to polar coordinates, where z points out of plane of orbit for each particle
            r_i,vr_i,phi_i,vphi_i,j_hat_i = xyz_to_polar(x_i,y_i,z_i,vx_i,vy_i,vz_i)
            r_f,vr_f,phi_f,vphi_f,j_hat_f = xyz_to_polar(x_f,y_f,z_f,vx_f,vy_f,vz_f)
        
        # r gets treated like a normal spatial coordinate (except cannot go <0; ignore for now?)
        r = interpolation_for_movie_pos(dt,delta_t, r_i,r_f,vr_i,vr_f, s_i,s_f,nm_i_f,nm_f_i,
            order=1,pos_interp_flag=pos_interp_flag,threshold_for_overshoot=1.e-10)
        # predict the 'final' phi, theta, knowing that these are actually periodic:
        # note for phi, need third-order: inexact matching can leave 'wedges' missing from orbits
        phi = interpolation_for_movie_pos(dt,delta_t, phi_i,phi_f,
            vphi_i,vphi_f, s_i,s_f,nm_i_f,nm_f_i, periodic=1,order=3,
            threshold_for_overshoot=0.25)
        x_t = r * np.cos(phi);  y_t = r * np.sin(phi);
        j_hat = interpolation_for_movie_jhat(dt,delta_t, j_hat_i,j_hat_f, s_i,s_f,nm_i_f,nm_f_i )

        #if(stars!=0): 
        #    pos_interp_flag = (r < 

        ## define the absolute 'phi' relative to an arbitrary (but fixed) 90-degree rotation of the angular momentum vector j
        x_jhat = 0.*j_hat
        x_jhat[:,0]=0.*j_hat[:,0]; x_jhat[:,1]=-j_hat[:,2]; x_jhat[:,2]=j_hat[:,1]
        jj = np.sqrt( x_jhat[:,0]*x_jhat[:,0] + x_jhat[:,1]*x_jhat[:,1] + x_jhat[:,2]*x_jhat[:,2] )
        for j in [0,1,2]: x_jhat[:,j] /= jj
        ## generate y-vector by cross-product a x b = c (gaurantees right-hand rule)
        y_jhat = cross( j_hat, x_jhat )
        jj = np.sqrt( y_jhat[:,0]*y_jhat[:,0] + y_jhat[:,1]*y_jhat[:,1] + y_jhat[:,2]*y_jhat[:,2] )
        for j in [0,1,2]: y_jhat[:,j] /= jj
        z_t = 0.0 * x_t
        if(polar_vs_cylindrical == False):
            z_t = interpolation_for_movie_pos(dt,delta_t, zc_i,zc_f,vzc_i,vzc_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)#,order=1)
        ## now project back into our lab coordinate system
        x = x_t*x_jhat[:,0] + y_t*y_jhat[:,0] + z_t*j_hat[:,0]
        y = x_t*x_jhat[:,1] + y_t*y_jhat[:,1] + z_t*j_hat[:,1]
        z = x_t*x_jhat[:,2] + y_t*y_jhat[:,2] + z_t*j_hat[:,2]

        ## use cartesian interpolation for some points
        x_c=interpolation_for_movie_pos(dt,delta_t, x_i,x_f,vx_i,vx_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)#,order=1)
        y_c=interpolation_for_movie_pos(dt,delta_t, y_i,y_f,vy_i,vy_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)#,order=1)
        z_c=interpolation_for_movie_pos(dt,delta_t, z_i,z_f,vz_i,vz_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)#,order=1)

        ## combine velocities into initial and final matched vectors
        Nmatched=x_i[s_i].size;
        Nunmatched_i=x_i[nm_i_f].size;
        Nunmatched_f=x_f[nm_f_i].size;
        u=np.zeros(Nmatched+Nunmatched_i+Nunmatched_f);
        z0=0.*u+1.e-10; dz=0.*u; vr=0.*u; vphi=0.*u; rr=0.*u; dr=0.*u; jj=0.*u; vz=0.*u; rvphi_i=r_i*vphi_i; rvphi_f=r_f*vphi_f;
        if (Nmatched>0):
            vr[0:Nmatched] = np.sqrt((vr_i[s_i]**2.+vr_f[s_f]**2.)/2.)
            vphi[0:Nmatched] = np.sqrt((rvphi_i[s_i]**2.+rvphi_f[s_f]**2.)/2.)
            rr[0:Nmatched] = np.sqrt((r_i[s_i]**2.+r_f[s_f]**2.)/2.)
            dr[0:Nmatched] = np.abs(r_i[s_i]-r_f[s_f]) / np.abs(r_i[s_i]+r_f[s_f])
            jj[0:Nmatched] = (j_hat_i[s_i,0]*j_hat_f[s_f,0] + j_hat_i[s_i,1]*j_hat_f[s_f,1] + j_hat_i[s_i,2]*j_hat_f[s_f,2]) / \
                np.sqrt((j_hat_i[s_i,0]*j_hat_i[s_i,0] + j_hat_i[s_i,1]*j_hat_i[s_i,1] + j_hat_i[s_i,2]*j_hat_i[s_i,2]) * \
                        (j_hat_f[s_f,0]*j_hat_f[s_f,0] + j_hat_f[s_f,1]*j_hat_f[s_f,1] + j_hat_f[s_f,2]*j_hat_f[s_f,2]))
            vz[0:Nmatched] = 0
            if(use_polar_interpolation==1):
                z0[0:Nmatched] = np.sqrt(zc_i[s_i]*zc_i[s_i] + zc_f[s_f]*zc_f[s_f])
                dz[0:Nmatched] = np.abs(zc_i[s_i] - zc_f[s_f])
            if(polar_vs_cylindrical == False):
                vz[0:Nmatched] = np.sqrt((vzc_i[s_i]**2.+vzc_f[s_f]**2.)/2.)
        if(Nunmatched_i>0):
            vr[Nmatched:Nmatched+Nunmatched_i] = vr_i[nm_i_f]
            vphi[Nmatched:Nmatched+Nunmatched_i] = rvphi_i[nm_i_f]
            rr[Nmatched:Nmatched+Nunmatched_i] = r_i[nm_i_f]
            rr[Nmatched:Nmatched+Nunmatched_i] = -1
            jj[Nmatched:Nmatched+Nunmatched_i] = 1
            vz[Nmatched:Nmatched+Nunmatched_i] = 1
        if(Nunmatched_f>0):
            vr[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = vr_f[nm_f_i]
            vphi[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = rvphi_f[nm_f_i]
            rr[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = r_f[nm_f_i]
            rr[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = -1
            jj[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = 1
            vz[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = 1
        ## attempt to guess what's 'rotation supported' and not, and use appropriate interpolation
        vr2 = vr*vr + vz*vz
        v2 = vr2 + vphi*vphi
        
        ## good on large scales
        no_rot_support = (z0>2.) | ( z0/rr > 0.6 ) | ( (dz/z0>2.)&(z0 > 1.) ) | ( rr == -1 ) | (jj < 0.8) | ( vr2 > v2/2. ) | ( rr > 12. ) # usually does ok, too much 'breathing' at late times
        no_rot_support_favorcartesian = no_rot_support
        
        ## good on small scales
        no_rot_support = (rr == -1) | (rr > 12.) | (set_dx_frame>60.) | ((dr/rr > 2.2)&(rr>2)&(set_dx_frame>30.)) ##| (jj < (-1.+2.*set_dx_frame/60.))
        favorcartesian = (np.random.RandomState(m_i.size+m_f.size).rand(rr.size) < np.exp(-(30./set_dx_frame)**2.))
        no_rot_support[favorcartesian] = no_rot_support_favorcartesian[favorcartesian]
        
        x[no_rot_support] = x_c[no_rot_support]
        y[no_rot_support] = y_c[no_rot_support]
        z[no_rot_support] = z_c[no_rot_support]
    else:
        x=interpolation_for_movie_pos(dt,delta_t, x_i,x_f,vx_i,vx_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)
        y=interpolation_for_movie_pos(dt,delta_t, y_i,y_f,vy_i,vy_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)
        z=interpolation_for_movie_pos(dt,delta_t, z_i,z_f,vz_i,vz_f, s_i,s_f,nm_i_f,nm_f_i,pos_interp_flag=pos_interp_flag)
                

    x = x + (cen_f[0]-cen_i[0]) * dt / delta_t + cen_i[0]
    y = y + (cen_f[1]-cen_i[1]) * dt / delta_t + cen_i[1]
    z = z + (cen_f[2]-cen_i[2]) * dt / delta_t + cen_i[2]

    h=interpolation_for_movie_scalar(dt,delta_t, h_i,h_f, s_i,s_f,nm_i_f,nm_f_i, hsml_rescale=1, pos_interp_flag=pos_interp_flag)
    m=interpolation_for_movie_scalar(dt,delta_t, m_i,m_f, s_i,s_f,nm_i_f,nm_f_i, mass_rescale=1, pos_interp_flag=pos_interp_flag)
    c=interpolation_for_movie_scalar(dt,delta_t, c_i,c_f, s_i,s_f,nm_i_f,nm_f_i)
    #c=np.log10(interpolation_for_movie_scalar(dt,delta_t, 10.**c_i,10.**c_f, s_i,s_f,nm_i_f,nm_f_i))
    zm=interpolation_for_movie_scalar(dt,delta_t, zm_i,zm_f, s_i,s_f,nm_i_f,nm_f_i)

    gc.collect()
    return x,y,z,m,h,c,zm;






def interpolation_for_movie_pos( dt, delta_t, u_i,u_f, vu_i,vu_f, \
    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i , periodic=0, order=3,
    pos_interp_flag=np.array([0.]), threshold_for_overshoot=0.25):

    #order = 1

    dt=np.array(dt); delta_t=np.array(delta_t); ## make sure casting is ok
    tau = dt/delta_t;
    Nmatched=u_i[sub_id_i].size;
    Nunmatched_i=u_i[nomatch_from_i_in_f].size;
    Nunmatched_f=u_f[nomatch_from_f_in_i].size;
    u=np.zeros(Nmatched+Nunmatched_i+Nunmatched_f);

    # first compute for objects in (i) with a match in (f)
    if (Nmatched>0):
        xi = np.copy(u_i[sub_id_i])  ; xf = np.copy(u_f[sub_id_f])  ;
        vi = np.copy(vu_i[sub_id_i]) ; vf = np.copy(vu_f[sub_id_f]) ;
        
        order_vec = 0.*xi + order
        order_tol = threshold_for_overshoot
        vdenom = 1. / (order_tol * (1.e-36 + np.abs(vf+vi)))
        xdenom = 1. / (order_tol * (1.e-36 + np.abs(xf+xi)))
        
        order_vec[(np.abs(vf-vi)*vdenom  > 1)] = 1
        
        if (periodic==1):
            ## u is periodic in 2pi, so estimate its final position and wrap coordinate
            ## guess where its going to be to get 'best estimate' of the number of revolutions
            #x_exp = xi + 0.5*(vf+vi)*delta_t # actually can get some errors when vf/i switch signs...
            x_exp = xi + vi*delta_t
            ## now pin the final location to the appropriate value
            xf = get_closest_periodic_value( x_exp, xf , xi)
        
        ## check for over/under-shooting, and if present, lower the reconstruction order
        #if(order==3):
        x0=xi; x1=vi*delta_t; x2=3.*(xf-xi)-(2.*vi+vf)*delta_t; x3=-2.*(xf-xi)+(vi+vf)*delta_t  
        # velocity extremum 
        t_x = -x2/(3.*x3)
        v_x=(x1-x2*t_x)/delta_t; 
        dv = np.minimum(np.abs(v_x-vi) , np.abs(v_x-vf))
        ok_tmp = (dv*vdenom > 1) & (t_x > 0.) & (t_x < 1.)
        order_vec[ok_tmp] = np.minimum(order_vec[ok_tmp],2 + np.zeros(order_vec[ok_tmp].size))
        # position extremum 1 (b/c quadratic solution)
        q=x2*x2 - 3.*x1*x3
        qq = t_x + np.sqrt(q)/(3.*x3)
        ok = (q > 0) & (qq > 0.) & (qq < 1.)
        xx=x0 + x1*qq + x2*qq*qq + x3*qq*qq*qq
        dx = np.minimum(np.abs(xx-xi) , np.abs(xx-xf))
        ok_tmp = (dx*xdenom > 1) & ok
        ok_tmp = (q > 0) & (qq > 0.) & (qq < 1.)
        order_vec[ok_tmp] = np.minimum(order_vec[ok_tmp],2 + np.zeros(order_vec[ok_tmp].size))
        # position extremum 2 (b/c quadratic solution)
        qq = t_x - np.sqrt(q)/(3.*x3)
        ok = (q > 0) & (qq > 0.) & (qq < 1.)
        xx=x0 + x1*qq + x2*qq*qq + x3*qq*qq*qq
        dx = np.minimum(np.abs(xx-xi) , np.abs(xx-xf))
        ok_tmp = (dx*xdenom > 1) & ok
        ok_tmp = (q > 0) & (qq > 0.) & (qq < 1.)
        order_vec[ok_tmp] = np.minimum(order_vec[ok_tmp],2 + np.zeros(order_vec[ok_tmp].size))

        #if(order==2):
        x0=xi; x2=(vf-vi)*delta_t/2.; x1=(xf-xi)-x2;
        t_x=-x1/(2.*x2)
        xx = x0 + (x1/2.)*t_x;
        dx = np.minimum(np.abs(xx-xi) , np.abs(xx-xf))
        ok_tmp = (dx*xdenom > 1) & (t_x > 0.) & (t_x < 1.) & (order_vec <= 2.)
        ok_tmp = (t_x > 0.) & (t_x < 1.) & (order_vec <= 2.)
        order_vec[ok_tmp] = np.minimum(order_vec[ok_tmp],1 + np.zeros(order_vec[ok_tmp].size))
        
        #order_vec = 0.*order_vec + 3.
        
        ## third-order interpolation: enables exact matching of x, v, but can 'overshoot'
        x0=xi; x1=vi*delta_t; x2=3.*(xf-xi)-(2.*vi+vf)*delta_t; x3=-2.*(xf-xi)+(vi+vf)*delta_t  
        u_o3 = x0 + x1*tau + x2*tau*tau + x3*tau*tau*tau
        ## second-order: minimizes absolute velocity difference (more stable, no exact velocity matching)
        x0=xi; x2=(vf-vi)*delta_t/2.; x1=(xf-xi)-x2;
        u_o2 = x0 + x1*tau + x2*tau*tau
        ## use linear interpolation below if for some reason this is unstable: least accurate, but most stable
        x0=xi; x1=(xf-xi)
        u_o1 = x0 + x1*tau
        
        #x0=xi; x1=vi*delta_t; x4=9.*(vf-vi)*delta_t/4.; 
        #x4=7.*(vf-vi)*delta_t/4.; 
        #x3=2.*(xi-(xf+x4))+(vf+vi)*delta_t; x2=0.5*((vf-vi)*delta_t-3.*x3-4.*x4);
        #u_o1 = x0 + x1*tau + x2*tau*tau + x3*tau*tau*tau + x4*tau*tau*tau*tau 
        
        #u_o1 = xi + (xf-xi)*tau + 0.5*((vf+vi)*delta_t+2.*(xi-xf))*np.sin(2.*np.pi*tau)/(2.*np.pi)
        
        u_tmp = u_o1
        order_vec=np.round(order_vec)
        ok = (order_vec == 3)
        u_tmp[ok] = u_o3[ok]
        ok = (order_vec == 2)
        u_tmp[ok] = u_o2[ok]
        #u[0:Nmatched] = u_tmp[0:Nmatched]
        
        #v0 = (xf-xi)/delta_t
        #dv = np.sqrt((vi-v0)**2 + (vf-v0)**2)
        #ok_tmp = (dv/(1.e-40 + np.abs(v0)) > 0.5) # too many at x_i, but does resolve central 'ring'
        #ok_tmp = (dv/(1.e-40 + np.abs(v0)) > 5.) # ring re-appears, potentially better?
        #ok_tmp = (dv/(1.e-40 + np.abs(v0)) > 1.6) # somewhere in-between, but not great
        #u_o1[ok_tmp] = xf[ok_tmp]
        #u_o1 = xf
        
        #dx = np.abs(xf-xi)
        #ok_tmp = (dx > 1.) # ok, good, but noisy
        #ok_tmp = (dx > 10.) # ring remains
        #ok_tmp = (dx > 3.) # ring remains
        
        if(pos_interp_flag.size > 1):
            ## smoothed step-function interpolation between snapshots
            width_ftau = 0.05
            width_ftau = 0.25
            q = 2.*np.pi/width_ftau
            f0 = 0.5-np.arctan(q*(0.-0.5))/np.pi
            f1 = 0.5-np.arctan(q*(1.-0.5))/np.pi
            f_tau = (0.5-np.arctan(q*(tau-0.5))/np.pi - f1) / (f0 - f1)
            sharp_interp = (xi + vi*dt) * f_tau + (xf + vf*(dt-delta_t)) * (1.-f_tau)
            #sharp_interp = (xi + 0*vi*dt) * f_tau + (xf + 0*vf*(dt-delta_t)) * (1.-f_tau)
            #u_tmp[pos_interp_flag] = sharp_interp[pos_interp_flag]
            #u_tmp[pos_interp_flag] = u_o3[pos_interp_flag]
            #if(tau < 0.5):
            #    u_tmp[pos_interp_flag] = xi[pos_interp_flag]
            #else:
            #    u_tmp[pos_interp_flag] = xf[pos_interp_flag]
            
        u[0:Nmatched] = u_tmp[0:Nmatched]
        
        
    ## now unmatched: first those in (i) with no match in (f)
    if(Nunmatched_i>0):
        u[Nmatched:Nmatched+Nunmatched_i] = \
            u_i[nomatch_from_i_in_f]+vu_i[nomatch_from_i_in_f]*dt;
	## now unmatched: now those in (f) with no match in (i)
    if(Nunmatched_f>0):
        u[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = \
            u_f[nomatch_from_f_in_i]+vu_f[nomatch_from_f_in_i]*(dt-delta_t);

    gc.collect()
    return u;


def interpolation_for_movie_jhat( dt, delta_t, j_i,j_f, \
    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i ):

    dt=np.array(dt); delta_t=np.array(delta_t); ## make sure casting is ok
    tau = dt/delta_t;
    u_i = np.zeros(j_i[:,0].size)
    u_f = np.zeros(j_f[:,0].size)
    Nmatched=u_i[sub_id_i].size;
    Nunmatched_i=u_i[nomatch_from_i_in_f].size;
    Nunmatched_f=u_f[nomatch_from_f_in_i].size;
    j_hat=np.zeros((Nmatched+Nunmatched_i+Nunmatched_f,3));

    # first compute for objects in (i) with a match in (f)
    if (Nmatched>0):
        ## initialize matrices
        pt_i = np.zeros(Nmatched); pt_f = 0.*pt_i; 
        p_i = np.zeros((Nmatched,3)); p_f = 0.*p_i; v_f = 0.*p_f
        for j in [0,1,2]:
            p_i[:,j] = j_i[sub_id_i,j]
            p_f[:,j] = j_f[sub_id_f,j]
            pt_i += p_i[:,j]*p_i[:,j]
            pt_f += p_f[:,j]*p_f[:,j]
        ## normalize vectors to unity
        pt_i = np.sqrt(pt_i); pt_f = np.sqrt(pt_f);
        for j in [0,1,2]:
            p_i[:,j] /= pt_i
            p_f[:,j] /= pt_f
            
        ## build a rotation matrix between the j_hat values at two different timesteps:
        cos_rot_ang = p_i[:,0]*p_f[:,0] + p_i[:,1]*p_f[:,1] + p_i[:,2]*p_f[:,2]
        angle = np.arccos( cos_rot_ang )

        angle_to_rotate = angle * tau 
        cos_ang = np.cos(angle_to_rotate)
        sin_ang = np.sin(angle_to_rotate)

        ## build perpendicular vectors and normalize
        k_rot = cross( p_i, p_f )
        kk=0.*pt_i;
        for j in [0,1,2]: kk += k_rot[:,j]*k_rot[:,j]
        kk = np.sqrt(kk)
        for j in [0,1,2]: k_rot[:,j] /= kk
        
        k_cross_p = cross( k_rot, p_i )
        kk=0.*pt_i;
        for j in [0,1,2]: kk += k_cross_p[:,j]*k_cross_p[:,j]
        kk = np.sqrt(kk)
        for j in [0,1,2]: k_cross_p[:,j] /= kk
        
        k_dot_p = k_rot[:,0]*p_i[:,0] + k_rot[:,1]*p_i[:,1] + k_rot[:,2]*p_i[:,2]
        for j in [0,1,2]:
            v_f[:,j] = p_i[:,j]*cos_ang + k_cross_p[:,j]*sin_ang + k_rot[:,j]*k_dot_p*(1.-cos_ang)

        for j in [0,1,2]:
            j_hat[0:Nmatched,j] = v_f[:,j]  ## now have the j_hat vector for this time

    ## now unmatched: first those in (i) with no match in (f)
    if(Nunmatched_i>0):
        for j in [0,1,2]:
            j_hat[Nmatched:Nmatched+Nunmatched_i,j] = j_i[nomatch_from_i_in_f,j]

	## now unmatched: now those in (f) with no match in (i)
    if(Nunmatched_f>0):
        for j in [0,1,2]:
            j_hat[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = j_f[nomatch_from_f_in_i]

    gc.collect()
    return j_hat;



def interpolation_for_movie_scalar(dt, delta_t, u_i,u_f, \
    sub_id_i,sub_id_f, nomatch_from_i_in_f,nomatch_from_f_in_i, \
    mass_rescale=0, hsml_rescale=0, pos_interp_flag=np.array([0.]) ):

    dt=np.array(dt); delta_t=np.array(delta_t); ## make sure casting is ok
    tau = dt/delta_t;
    Nmatched=u_i[sub_id_i].size;
    Nunmatched_i=u_i[nomatch_from_i_in_f].size;
    Nunmatched_f=u_f[nomatch_from_f_in_i].size;
    u=np.zeros(Nmatched+Nunmatched_i+Nunmatched_f);
    
    # first compute for objects in (i) with a match in (f)
    if (Nmatched>0):
        #u[0:Nmatched] = u_i[sub_id_i] + (u_f[sub_id_f]-u_i[sub_id_i])*tau
        u[0:Nmatched] = u_i[sub_id_i] + (u_f[sub_id_f]-u_i[sub_id_i])*tau
        
        u0 = np.maximum(u_i[sub_id_i],1.e-10)
        u1 = np.maximum(u_f[sub_id_f],1.e-10)
        u_tmp = u0 * np.exp(tau * np.log(u1/u0))
        u[0:Nmatched] = u_tmp 

        if(hsml_rescale==1):
            u0 = 1.*u_i[sub_id_i]; u1 = 1.*u_f[sub_id_f];
            u_tmp = u0*(1.-tau) + u1*tau
            if(pos_interp_flag.size > 1):
                u0=u0[pos_interp_flag]; u1=u1[pos_interp_flag];
                u_tmp[pos_interp_flag] = np.sqrt(u0*u0*(1.-tau) + u1*u1*tau);
            u[0:Nmatched] = u_tmp
        if(mass_rescale==1):
            u_tmp = u_i[sub_id_i] + (u_f[sub_id_f]-u_i[sub_id_i])*tau
            if(pos_interp_flag.size > 1):
                depth = 0.1 # fraction to diminish, which sets width of dip
                #u_tmp_1 = 1.*u_tmp; u_tmp *= 0.;
                u_tmp_1 = 0.*u_tmp; 
                #u_tmp[pos_interp_flag] = u_tmp_1[pos_interp_flag]
                #u_tmp[pos_interp_flag] /= (1. + tau*(1.-tau) * 4.*(1.-depth)/depth)
            u[0:Nmatched] = u_tmp
            
    ## now unmatched: first those in (i) with no match in (f)
    if(Nunmatched_i>0):
        if(mass_rescale==1):
            u_t = 1.*u_i[nomatch_from_i_in_f] * ((1.-tau)**(2.))
        elif(hsml_rescale==1):
            u_t = u_i[nomatch_from_i_in_f] / ((1.-tau)**2.+1.0e-5)
        else:
            u_t = u_i[nomatch_from_i_in_f]
        u[Nmatched:Nmatched+Nunmatched_i] = u_t;
	## now unmatched: now those in (f) with no match in (i)
    if(Nunmatched_f>0):
        if(mass_rescale==1):
            #u_t = u_f[nomatch_from_f_in_i] * (tau**(2.))
            u_t = 1.*u_f[nomatch_from_f_in_i] * (tau**(3.))
        elif(hsml_rescale==1):
            u_t = u_f[nomatch_from_f_in_i] / ((1.-(1.-tau))**2.+1.0e-5)
        else:
            u_t = u_f[nomatch_from_f_in_i]
        u[Nmatched+Nunmatched_i:Nmatched+Nunmatched_i+Nunmatched_f] = u_t;

    gc.collect()
    return u;
