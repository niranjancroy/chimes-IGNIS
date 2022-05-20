import numpy as np
import matplotlib
#matplotlib.use('Agg') ## this calls matplotlib without a X-windows GUI
from subprocess import call
import matplotlib.colors
import matplotlib.cm
import matplotlib.pyplot as plt
import math
import pfh_utils as util
import gadget
import gadget_lib.load_stellar_hsml as starhsml
import visualization.colors as viscolors
import visualization.get_attenuated_stellar_luminosities as getlum
import visualization.make_threeband_image as makethreepic
import visualization.contour_makepic as cmakepic
import visualization.raytrace_projection as rayproj
import h5py


def checklen(x):
    return len(np.array(x,ndmin=1));

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (np.abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (np.abs(input)<=xmax);

def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;
	
def get_gas_xray_luminosity(ppp):
    brems = gadget.gas_xray_brems( \
        ppp['m'], ppp['u'], ppp['rho'], ppp['ne'], ppp['nh'] );
    return brems; ## can also add metal line cooling luminosities later

##
## handy function to load everything we need from snapshots, in particular to 
##   concatenate the various star types if we're using non-cosmological snapshots
##
def load_snapshot_brick(ptypes, snapdir, snapnum, h0=1, cosmological=0, \
        skip_bh=0, do_xray=0, four_char=0, use_rundir=0, full_gas_set=0 , filename_prefix=''):
    have=0; have_h_stars=0;
    ppp_head=gadget.readsnap(snapdir,snapnum,0,h0=h0,cosmological=cosmological,skip_bh=skip_bh,header_only=1)
    time=ppp_head['time']
    if (full_gas_set==1):
        ## here it will just return the entire gas 'brick'
        ppp=gadget.readsnap(snapdir,snapnum,0,h0=h0,cosmological=cosmological,skip_bh=skip_bh);
        m=ppp['m']; p=ppp['p']; x=p[:,0]; y=p[:,1]; z=p[:,2]; zm=ppp['z']; 
        if(checklen(zm.shape)>1): zm=zm[:,0]; lx=get_gas_xray_luminosity(ppp) / (3.9e33);
        sfr0=ppp['SFR']; sfr=1.*sfr0; 
        ok=sfr0[(sfr0 > 0.)]
        sfr = 1.e-10 * 0.25*(1./(1.+(ppp['rho'])**(-0.5)))
        if(ok.size > 0):
            p0=np.min(ppp['rho'][(sfr0 > 0.)]); 
            lo=(sfr0 <= 0.); sfr[lo]=0.25*(1./(1.+(ppp['rho'][lo]/p0)**(-0.5)))*np.min(sfr0[(sfr0>0.)]);
        return ppp['u'], ppp['rho'], ppp['h'], ppp['nh'], ppp['ne'], sfr, lx, zm, m, x, y, z, time;
        
    for ptype in ptypes:
        ppp=gadget.readsnap(snapdir,snapnum,ptype,h0=h0,cosmological=cosmological,skip_bh=skip_bh);
        if(ppp['k']==1):
            n=checklen(ppp['m']);
            if(n>1):
                m=ppp['m']; p=ppp['p']; x=p[:,0]; y=p[:,1]; z=p[:,2]; 
                if (ptype==0): 
                    cc=gadget.gas_temperature(ppp['u'],ppp['ne']);
                    if (do_xray==1): m=get_gas_xray_luminosity(ppp) / (3.9e33); ## bremstrahhlung (in solar)
                    h=ppp['h']; zm=ppp['z']; 
                    if(checklen(zm.shape)>1): zm=zm[:,0]
                if (ptype==4):
                    cc=gadget.get_stellar_ages(ppp,ppp_head,cosmological=cosmological);
                    zm=ppp['z']; 
                    if(checklen(zm.shape)>1): zm=zm[:,0]
                if (ptype==2): ## need to assign ages and metallicities
                    cc=np.random.rand(n)*(ppp_head['time']+4.0);
                    zm=(np.random.rand(n)*(1.0-0.1)+0.1) * 0.02;
                if (ptype==3): ## need to assign ages and metallicities
                    cc=np.random.rand(n)*(ppp_head['time']+12.0);
                    zm=(np.random.rand(n)*(0.3-0.03)+0.03) * 0.02;
                hstars_should_concat=0;
                if((ptype>0) & (have_h_stars==0)):
                    h=starhsml.load_allstars_hsml(snapdir,snapnum,cosmo=cosmological, \
                        use_rundir=use_rundir,four_char=four_char,use_h0=h0,filename_prefix=filename_prefix);
                    have_h_stars=1;
                    hstars_should_concat=1;
                if (have==1):
                    m_all=np.concatenate((m_all,m)); x_all=np.concatenate((x_all,x)); 
                    y_all=np.concatenate((y_all,y)); z_all=np.concatenate((z_all,z)); 
                    c_all=np.concatenate((c_all,cc)); zm_all=np.concatenate((zm_all,zm)); 
                    if(hstars_should_concat==1): h_all=np.concatenate((h_all,h)); 
                else:
                    m_all=m; x_all=x; y_all=y; z_all=z; zm_all=zm; c_all=cc; h_all=h; 
                    have=1;
    if (have==1):
        return m_all, x_all, y_all, z_all, c_all, h_all, zm_all, time;
    return 0,0,0,0,0,0,0,0;
    

def overlay_scale_label(xr,yr,figure_axis,c0='w'):
    ddx=0.5*(xr[1]-xr[0]);
    if (ddx>0.00002): dx=0.00001; dxl='0.01 pc';
    if (ddx>0.0002): dx=0.0001; dxl='0.1 pc';
    if (ddx>0.002): dx=0.001; dxl='1 pc';
    if (ddx>0.02): dx=0.01; dxl='10 pc';
    if (ddx>0.2): dx=0.1; dxl='100 pc';
    if (ddx>2.): dx=1.; dxl='1 kpc';
    if (ddx>20.): dx=10.; dxl='10 kpc';
    if (ddx>200.): dx=100.; dxl='100 kpc';
    if (ddx>2000.): dx=1000.; dxl='1 Mpc';
    if (ddx>20000.): dx=10000.; dxl='10 Mpc';
    
    xlen=(xr[1]-xr[0])
    ylen=(yr[1]-yr[0])
    xoff = (0.25+0.02)*ddx / xlen
    yoff = 0.025*(yr[1]-yr[0]) / ylen
    xr_new = np.array([xoff-dx/xlen*0.5,xoff+dx/xlen*0.5])
    yr_new = np.array([yoff,yoff])
    
    plt.text(xoff,1.75*yoff,dxl,color=c0,\
        horizontalalignment='center',verticalalignment='baseline',\
        transform=figure_axis.transAxes,fontsize=3)
    figure_axis.autoscale(enable=False,axis='both',tight=True)
    figure_axis.plot(xr_new,yr_new,color=c0,linewidth=0.7,transform=figure_axis.transAxes)

    
def overlay_time_label(time, figure_axis, cosmo=0, c0='w', n_sig_add=0, tunit_suffix='Gyr'):
    if(cosmo==0): 
        time_to_use = time ## absolute time in Gyr
        prefix = ''
        suffix = tunit
    else: ## 
        time_to_use = 1./time - 1. ## redshift
        prefix = 'z='
        suffix = ''
        
    n_sig = 2
    if time_to_use >= 10: n_sig += 0
    if time_to_use < 1.: n_sig += 1
    t_str = round_to_n( time_to_use, n_sig+n_sig_add )
    label_str = prefix+t_str+suffix
    xoff=0.03; yoff=1.-0.025;
    plt.text(xoff,yoff,label_str,color=c0,\
        horizontalalignment='left',verticalalignment='top',\
        transform=figure_axis.transAxes,fontsize=3)


def round_to_n(x, n):
    ''' Utility function used to round labels to significant figures for display purposes  from: http://mail.python.org/pipermail/tutor/2004-July/030324.html'''
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")

    # show everything as floats (preference; can switch using code below to showing eN instead
    format = "%." +str(n-1) +"f"
    as_string=format % x
    return as_string

    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    if as_string[-3:] in ['+00', '+01', '+02', '+03','-01', '-02', '-03']:
        #then number is 'small', show this as a float
        format = "%." +str(n-1) +"f"
        as_string=format % x
    return as_string


def coordinates_rotate(x_all, y_all, z_all, theta, phi, coordinates_cylindrical=0):
    ## set viewing angle for image
    ## project into plane defined by vectors perpendicular to the vector of this angle
    x = np.cos(phi)*x_all + np.sin(phi)*y_all + 0.*z_all
    y = -np.cos(theta)*np.sin(phi)*x_all + np.cos(theta)*np.cos(phi)*y_all + np.sin(theta)*z_all
    z =  np.sin(theta)*np.sin(phi)*x_all - np.sin(theta)*np.cos(phi)*y_all + np.cos(theta)*z_all
    if (coordinates_cylindrical==1):
        x=x/abs(x)*np.sqrt(x*x+z*z); y=y; z=z
    return x,y,z
  
 
#
# given a set of positions, project them to the coordinates they 'would be at'
#   with respect to a given camera (for ray-tracing purposes)
#
def coordinates_project_to_camera(x, y, z, \
        camera_pos=[0.,0.,0.], camera_dir=[0.,0.,1.], \
        screen_distance=1.0 ):
    camera_pos=np.array(camera_pos,dtype='d');
    camera_dir=np.array(camera_dir,dtype='d');
    ## recenter at camera location
    x-=camera_pos[0]; y-=camera_pos[1]; z-=camera_pos[2];
    ## line-of-sight distance along camera viewing angle
    znew = x*camera_dir[0] + y*camera_dir[1] + z*camera_dir[2]; 
    ## get perpendicular axes to determine position in 'x-y' space:
    cperp_x,cperp_y = util.return_perp_vectors(camera_dir);
    xnew = x*cperp_x[0] + y*cperp_x[1] + z*cperp_x[2];
    ynew = x*cperp_y[0] + y*cperp_y[1] + z*cperp_y[2];

    ## project to 'screen' at fixed distance along z-axis :: 
    znew[znew==0.]=-1.0e-10; ## just to prevent nans
    ## need correction factors for other quantities (like hsml, for example)
    r_corr = screen_distance/znew;
    xnew *= r_corr; ## put positions 'at the screen'
    ynew *= r_corr; ## put positions 'at the screen'
    znew = -znew; ## so values in the direction of camera are *negative* (for later sorting)

    return xnew, ynew, znew, r_corr


def return_rotation_matrix(theta):
   tx,ty,tz = theta
   Rx = np.array([[1,0,0], [0, np.cos(tx), -np.sin(tx)], [0, np.sin(tx), np.cos(tx)]])
   Ry = np.array([[np.cos(ty), 0, -np.sin(ty)], [0, 1, 0], [np.sin(ty), 0, np.cos(ty)]])
   Rz = np.array([[np.cos(tz), -np.sin(tz), 0], [np.sin(tz), np.cos(tz), 0], [0,0,1]])
   return np.dot(Rx, np.dot(Ry, Rz))


def get_snap_cen(sdir):
    cen=[0.,0.,0.]
    if(sdir=='m09_ref11_adapt'): cen=[ 4710.58752192,  3425.17425247,  3176.85378332]
    if(sdir=='m09_ref11_adapt_metaldiffmhd_hiThold'): cen=[ 4709.3307902,   3425.09992071,  3178.07732329]
    if(sdir=='m09_ref11_fixedsoft'): cen=[ 4707.89059383, 3424.72304372, 3178.41580743]
    if(sdir=='sf_model_tests_default'): cen=[ 41828.52813216,  44147.23738634,  46238.4724623 ]
    if(sdir=='m12i_ref11'): cen=[ 41718.11654436,  44229.69343256,  46423.92902392]
    if(sdir=='m12i_ref11_tsph'): cen=[ 41817.80042241,  44157.58171149,  46286.9059771 ]
    if(sdir=='m10v_ref8'):  cen=[ 3403.30284402,  3298.51326896,  3529.32402938]
    if(sdir=='m10v_ref9'):  cen=[ 3397.61109927,  3290.66750332,  3524.44558671]
    if(sdir=='m10v_ref10'): cen=[ 3418.28945977,  3291.88373753,  3566.75524528]
    if(sdir=='m10v_ref11_adapt'): cen=[ 3427.06370651,  3271.65994972,  3523.04310327]
    if(sdir=='m10v_ref11_adapt_metaldiffmhd'): cen=[ 3426.16529536,  3271.85564295,  3522.90864404]
    if(sdir=='m10v_ref11'): cen=[ 3427.06370651,  3271.65994972,  3523.04310327]
    if(sdir=='m10q_ref8'):  cen=[ 3335.72293591,  2974.47073640,  3350.93201216]
    if(sdir=='m10q_ref9'):  cen=[ 3329.01555615,  2968.60596968,  3343.08745896]
    if(sdir=='m10q_ref10'): cen=[ 3329.05788868,  2985.49113833,  3311.77384541]
    if(sdir=='m10q_ref11'): cen=[ 3313.02295174,  2945.93301173,  3394.93218315]
    if(sdir=='m10q_ref11_fixed'): cen=[ 3316.39482422, 2943.0249495,  3393.59247881]
    if(sdir=='m10q_ref11_adapt_alldiffmhd'): cen=[ 3303.94783806,  2947.05529221,  3397.22271046]
    if(sdir=='m10q_ref11_adapt_hiThold'): cen=[ 3313.13005136,  2943.92293291,  3397.17398351]
    if(sdir=='m12i_ref11'): cen=[ 41718.11654436,  44229.69343256,  46423.92902392]
    if(sdir=='m12i_ref12'): cen=[ 41828.29846686,  44131.82912578,  46298.30386540]
    if(sdir=='m12i_ref13'): cen=[ 41875.75676213,  44122.38220725,  46257.48356735]
    if(sdir=='m12i_ref12_adapt_metaldiff2'): cen=[41830.50437144,  44239.39314392,  46234.6609263]
    if(sdir=='m12i_ref11_noRad_tsph'): cen=[ 41818.87257281,  44160.69550523,  46289.22666334]
    if(sdir=='m12f_ref12_fb-sym'): cen=[ 38748.69056642,  47644.00698554,  46763.80678027]
    if(sdir=='m12f_ref12_fb-iso'): cen=[ 38680.5433167,   47707.30194236,  46818.98373489]
    if(sdir=='m12b_ref12_fb-sym'): cen=[ 39311.00446494, 41979.01997792, 39152.44726635]
    if(sdir=='m12b_ref12_fb-iso'): cen=[ 39337.51460655,  41756.55378332,  39179.1354953 ]
    if(sdir=='m12c_ref12_fb-sym'): cen=[ 36008.06717868,  49152.75510195,  46821.18188375]
    if(sdir=='m12c_ref12_fb-iso'): cen=[ 35916.69167078,  49290.06782276,  46877.39851561]
    if(sdir=='m12m_ref12_fb-sym'): cen=[ 37561.37411999,  41583.64931681,  46769.30807611]
    if(sdir=='m12m_ref12_fb-iso'): cen=[ 37553.97725373,  41653.68221635,  46756.85800067]
    if(sdir=='m12f_ref13_fb-sym'): cen=[ 38744.24350156,  47701.65789227,  46790.76338099]
    if(sdir=='m11v_ref13_adapt'): cen=[ 42824.85443495,  44154.3228809,   43528.57110144]
    if(sdir=='m11q_ref13_adapt'): cen=[ 42242.33758432,  45740.82325377,  40512.26357658]
    if(sdir=='m12i_ref12_fire1'): cen=[ 41925.08412405,  44677.04569901,  46068.99276398]
    if(sdir=='m09_ref11_fire1'): cen=[ 4705.95809371,  3424.56268181,  3178.0762833 ]
    if(sdir=='m10q_ref11_fire1'): cen=[ 3207.62480088,  2976.60457247,  3286.54402602]
    if(sdir=='m10v_ref10_fire1'): cen=[ 4494.65998758,  3167.99864101,  2686.66612086]
    if(sdir=='m11v_ref12_fire1'): cen=[ 43229.3474205,   44158.01019037,  43678.13068558]
    if(sdir=='m11q_ref13_fire1'): cen=[ 42051.52966859,  45901.56828362,  39757.28554621]
    if(sdir=='m12i_ref12_fb_iso'): cen=[ 41843.52124705,  44171.67922194,  46253.53395134]
    if(sdir=='m12i_ref11_fb-aniso'): cen=[ 41710.99239422,  44338.27894807,  46346.86432687]
    if(sdir=='m12i_ref11_fb_iso'): cen=[ 41797.51755388,  44322.46511334,  46321.86212181]
    if(sdir=='m12i_ref12_fb_aniso'): cen=[ 41867.35752106,  44164.56941237,  46238.50485002]
    if(sdir=='m10q_ref10_fixedsoft_20'): cen=[ 3330.1619776,   2981.91197077,  3312.14271496]
    if(sdir=='m10q_ref10_fixedsoft_dm'): cen=[ 3318.74676144,  2982.3481015,   3316.68912189]    
    if(sdir=='m10q_ref11_fixedsoft'): cen=[ 3316.39482422,  2943.0249495,   3393.59247881]
    if(sdir=='m10q_ref11_dm_fixedsoft'): cen=[ 3302.74423667,  2947.34243833,  3400.34404091]
    if(sdir=='m10v_ref10_fixedsoft_dm'): cen=[ 3407.94699333,  3299.88783958, 3559.30869286]
    if(sdir=='m10v_ref11_fixedsoft'): cen=[ 3426.82879833,  3271.68589105,  3522.73985386]
    if(sdir=='m10v_ref11_dm_fixedsoft'): cen=[ 3415.57177299,  3277.97178013,  3517.26401981]
    if(sdir=='m12i_ref12_dm'): cen=[ 41824.99971528,  44152.94674268,  46280.59357175]
    if(sdir=='m12i_ref13_dm'): cen=[ 41820.51225595,  44153.30746753,  46272.85728992]
    if(sdir=='m12f_ref11_fb-sym'): cen=[ 38616.98367136,  47776.05817846,  46847.70543562]
    if(sdir=='m12f_ref12_dm'): cen=[ 38723.17211486,  47653.0094322,   46787.35209018]
    if(sdir=='m12f_ref12_dm_rad6'): cen=[ 38724.56059854,  47658.2378309,   46786.7506989 ]
    if(sdir=='m12f_ref13_dm'): cen=[ 38715.93421233,  47652.45341032,  46782.89188415]
    if(sdir=='m11q_ref13_fixedsoft'): cen=[ 42259.63436309,  45737.67365735,  40502.37933148]
    if(sdir=='m11v_ref13_fixedsoft'): cen=[ 42832.21307527,  44182.01152496,  43521.16148749]
    if(sdir=='m12f_ref11_dm_fixedsoft'): cen=[ 38732.23485065,  47642.62676188,  46775.60105084]
    if(sdir=='m12i_ref11_dm_smallsoft'): cen=[ 41827.37934501,  44152.90244138,  46274.98437919]
    if(sdir=='m12i_ref11_dm_small2soft'): cen=[ 41827.73049651,  44151.33876179,  46274.76064483]
    if(sdir=='m12i_ref11_dm'): cen=[ 41827.3526075,   44151.2047719,   46274.29765069]
    if(sdir=='m12i_ref11_dm_bigsoft'): cen=[ 41828.11838077,  44150.94670074,  46274.22638959]
    if(sdir=='m12i_ref11_dm_comoving600pc'): cen=[ 41826.84372491,  44152.31257079,  46274.1704996 ]
    if(sdir=='m12i_ref11_dm_ags'): cen=[ 41828.41312479,  44152.51080231,  46272.72380872]
    if(sdir=='m12i_ref13_dm_res-20pc'): cen=np.array([ 41820.82506815,  44152.98201607, 46272.80317769])
    if(sdir=='m12i_ref13_dm_res-80pc'): cen=np.array([ 41820.71264973,  44153.48068267,  46272.7722742 ])
    if(sdir=='m12i_ref13_dm_res-160pc'): cen=np.array([ 41820.83440886,  44153.3748299,   46272.86645861])
    if(sdir=='m10v_ref8_dm'): cen=[ 3399.08016416,  3303.44614895,  3522.91628842]
    if(sdir=='m10v_ref8_fixedsoft_dm'): cen=[ 3396.16711927,  3301.70932492,  3525.64239658]
    if(sdir=='m10v_ref9_dm'): cen=[ 3388.88283372,  3297.18020894,  3517.89061786]
    if(sdir=='m10v_ref9_fixedsoft_dm'): cen=[ 3388.85220596,  3297.1886898,   3517.74826294]
    if(sdir=='m10v_ref10_dm'): cen=[ 3408.40547067,  3299.4954469,   3558.31712581]
    if(sdir=='m10v_ref10_fixedsoft_dm'): cen=[ 3407.94699333,  3299.88783958,  3559.30869286]
    if(sdir=='m10v_ref11_dm_adapt'): cen=[ 3415.89241673,  3277.89113139,  3516.8104357 ]
    if(sdir=='m12i_ref13_fb_aniso'): cen=[ 41875.12718359,  44112.15174466,  46273.36882759]
    if(sdir=='m12i_ref13_MHD_CV_TD'): cen=[41815.491,44136.821,46267.662]
    return cen



def construct_spin_movie_single_sim(
    smaster = '/scratch/01799/phopkins/m12i_hybrid_test/m12_awetzel/m12i/fb-sym/', \
    omaster = '/work/01799/phopkins/images/', \
    sdir = 'm12i_ref12', \
    snum = 600, \
    key = 'star', \
    nf_base = 100, \
    add_key = '', \
    x00_in = 2.0, \
    x00_med = 30.0, \
    frac_to_rotate_initial_spin = 0.76, \
    pixels = 800., \
    theta_initial = 11., \
    theta_median = 109., \
    phi_initial = 0., \
    phi_median = -60., \
    frame_min = 0, \
    frame_max = -1,\
    scattered_fraction = 0.05,\
    h_rescale_factor=1.,\
    h_max=0):
    
    cen=get_snap_cen(sdir); print('cen = ',cen)
    if(cen[0]==0.):
        cen=gadget.calculate_zoom_center(smaster+'/'+sdir+'/output/',snum)
        print('CALCULATED CENTER == ',cen)
    
    P=gadget.readsnap(smaster+'/'+sdir+'/output/',snum,0,cosmological=1)
    x=P['p'][:,0]-cen[0]; y=P['p'][:,1]-cen[1]; z=P['p'][:,2]-cen[2]; r=np.sqrt(x*x+y*y+z*z)
    ok=np.where(r < 15.)
    x=x[ok]; y=y[ok]; z=z[ok]; r=r[ok]; m=P['m'][ok]; 
    vx=P['v'][ok,0]; vy=P['v'][ok,1]; vz=P['v'][ok,2];
    vx-=np.median(vx); vy-=np.median(vy); vz-=np.median(vz);
    jx=vy*z-vz*y; jy=vz*x-vx*z; jz=vx*y-vy*x;
    jxm=np.sum(m*jx); jym=np.sum(m*jy); jzm=np.sum(m*jz);
    jmm=np.sqrt(jxm*jxm+jym*jym+jzm*jzm)
    jxm/=jmm; jym/=jmm; jzm/=jmm;
    z_vec = np.array([jxm,jym,jzm])
    
    # now set up the spin
    nf1 = nf_base
    nf2 = nf_base
    nf2x = np.int(np.round(nf_base/1.5))
    nf3 = nf2
    frac1 = frac_to_rotate_initial_spin
    t0=theta_initial; t1=theta_median; p0=phi_initial; p1=phi_median;
    nf2 = np.int(np.round(nf_base*1.5))
    nf2x = np.int(np.round(nf_base*1.0))
    nf3 = nf_base
    
    n = 2.*nf1
    dt = frac1*(t1-t0)/n; dt0=dt;
    dp = frac1*(p1-p0)/n; dp0=dp;
    dlnx = 2.5/n; dlnx0=dlnx;
    theta_grid = t0 + dt*np.arange(0,n)
    phi_grid = p0 + dp*np.arange(0,n)
    x00_grid = x00_med * np.exp(dlnx*(n-np.arange(0,n)))

    n = nf2
    dt = 1.*(1.-frac1)*(t1-t0)/n;
    dp = (1.-frac1)*(p1-p0)/n;
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,x00_in,dlnx0,n)

    n = nf2x
    dt = 1.*(1.-frac1)*(t1-t0)/n
    dp = 2.*(1.-frac1)*(p1-p0)/n
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,0.5*x00_med,dlnx0,n)

    n = nf1
    dt = (frac1-1.*(1.-frac1))*frac1*(t1-t0)/n
    dp = (frac1-2.*(1.-frac1))*(p1-p0)/n
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,0.8*x00_med,dlnx0,n)

    n = nf3
    dt = 0.33*frac1*(t1-t0)/n
    dp = 0.5*180.0/n
    xf = 1.2

    n *= 2.5
    dt *= 1.5
    dp *= 1.5
    xf *= 5.

    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,xf*x00_med,dlnx0,n)

    print('ESTIMATED_NUMBER_OF_FRAMES == ',x00_grid.size)

    if(frame_max<0): frame_max=x00_grid.size-1;
    if(frame_min<0): frame_min=0;
    i_loop_list = np.arange(frame_min,frame_max+1,1)
    for i in i_loop_list:
        x00=x00_grid[i]; theta_0=theta_grid[i]; phi_0=phi_grid[i]; 
        addlayer='';set_added_layer_alpha=0.0;lightk=0;dustgas=1;ctable='heat_purple'
        threecolor=1; invert_colors=0;
        if(key=='star'):
            dynr = 2.e3; maxd = 1.1e-2 * (x00 / 30.)**(-1.5); lightk=0; 
            dustgas=2.2; 
            
            if (add_key=='nodust'): 
                dustgas = 1.e-5
            if (add_key=='CO'): 
                addlayer='CO'; ctable='heat_blue'; set_added_layer_alpha=0.4;
            if (add_key=='Halpha'): 
                addlayer='Halpha'; ctable='heat_green'; set_added_layer_alpha=0.3;
            if (add_key=='SFR'):
                addlayer='SFR'; ctable='heat_purple'; set_added_layer_alpha=0.3;
            if (add_key=='Xray'):
                addlayer='Xray'; ctable='heat_redyellow'; set_added_layer_alpha=0.4;
            if (add_key=='Zmetal'):
                addlayer='Zmetal'; ctable='heat_orange'; set_added_layer_alpha=0.8;
        if(key=='gas'):
            dynr = 1.e3; maxd = 6.0e-3 * (x00 / 30.)**(-1.5);
            threecolor=0; maxd*=10.0; dynr*=3.; lightk=0; invert_colors=0;
            #threecolor=1; lightk=1; ## not good for edge-on projections inside disk

        added_label='_full'
        if(add_key!=''): added_label='_'+add_key
        fname=omaster+'/'+sdir+'/'+sdir+'_'+key+added_label+'_f'+frame_ext(i)

	#theta = [phi_0, 0 , theta_0]
        #z_vec = np.dot(return_rotation_matrix(theta),z_vec_0)
	#z_vec = z_vec_0
        x_vec,y_vec = util.return_perp_vectors(z_vec)

        out0,out1 = image_maker(sdir,snum,center=cen,
            snapdir_master=smaster,outdir_master=omaster,pixels=pixels,
            do_with_colors=1,threecolor=threecolor, project_to_camera=1,cosmo=1,
            xrange=np.array([-1.,1.])*x00,yrange=np.array([-1.,1.])*x00,
            dynrange=dynr,maxden=maxd, show_gasstarxray=key,include_lighting=lightk,
            theta=theta_0,phi=phi_0,filename_set_manually=fname,set_added_layer_ctable=ctable,
            add_gas_layer_to_image=addlayer,dust_to_gas_ratio_rescale=dustgas,set_added_layer_alpha=set_added_layer_alpha,
            invert_colors=invert_colors,
            show_scale_label=0,show_time_label=0,spin_movie=1,
            scattered_fraction=scattered_fraction, BAND_IDS = [9,10,11], 
            h_rescale_factor=h_rescale_factor,h_max=h_max,
            projection_vector_x=x_vec,projection_vector_y=y_vec,projection_vector_z=z_vec);        






def image_maker( sdir, snapnum, \
    snapdir_master='/n/scratch2/hernquist_lab/phopkins/sbw_tests/', 
    outdir_master='/n/scratch2/hernquist_lab/phopkins/images/', 
    snapshot_subdirectory='output',
    theta=0., phi=0., dynrange=1.0e5, maxden_rescale=1., maxden=0.,
    show_gasstarxray = 'gas', #	or 'star' or 'xray'
    add_gas_layer_to_image='', set_added_layer_alpha=0.3, set_added_layer_ctable='heat_purple', 
    add_extension='', show_time_label=1, show_scale_label=1, 
    filename_set_manually='',
    do_with_colors=1, log_temp_wt=1, include_lighting=1, 
    set_percent_maxden=0, set_percent_minden=0, 
    center_on_com=0, center_on_bh=0, use_h0=1, cosmo=0, 
    center=[0., 0., 0.], coordinates_cylindrical=0, 
    offset=[0.,0.,0.],
    pixels=720,xrange=[-1.,1.],yrange=0,zrange=0,set_temp_max=0,set_temp_min=0,
    set_size_inches=1.,axis_enabled=False,
    threecolor=1, nasa_colors=1, sdss_colors=0, use_old_extinction_routine=0, 
    project_to_camera=1, camera_opening_angle=45.0,
    center_is_camera_position=0, camera_direction=[0.,0.,-1.], 
    gas_map_temperature_cuts=[1.0e4, 1.0e6], 
    input_data_is_sent_directly=0,
    m_all=0,x_all=0,y_all=0,z_all=0,c_all=0,h_all=0,zm_all=0,
    gas_u=0,gas_rho=0,gas_hsml=0,gas_numh=0,gas_nume=0,gas_metallicity=0,
    gas_mass=0,gas_x=0,gas_y=0,gas_z=0,time=0,
    invert_colors=0,spin_movie=0,scattered_fraction=0.01,
    min_stellar_age=0,h_rescale_factor=1.,h_max=0,
    angle_to_rotate_image=0.,
    projection_vector_x=[1,0,0],projection_vector_y=[0,1,0],projection_vector_z=[0,0,1],
    BAND_IDS=[9,10,11], ## ugr composite
    set_f_saturated=0.0001,
    dust_to_gas_ratio_rescale=1.0):
    

	## define some of the variables to be used
    ss=snap_ext(snapnum,four_char=1);
    tt=snap_ext(np.around(theta).astype(int));
    theta *= math.pi/180.; phi *= math.pi/180.; # to radians
    nameroot = sdir+'_s'+ss+'_t'+tt;
    outputdir = outdir_master+sdir;
    call(["mkdir",outputdir]);
    outputdir+='/'; 
    snapdir=snapdir_master+sdir+'/'+snapshot_subdirectory+'/';
    suff='_'+show_gasstarxray;
    if (threecolor==1):
        if (sdss_colors==1): suff+='_S3c'
        if (nasa_colors==1): suff+='_N3c'
    suff+=add_extension;
    fname_base=outputdir+nameroot+suff
    do_xray=0; do_stars=0; 
    if((show_gasstarxray=='xr') or (show_gasstarxray=='xray')): do_xray=1; do_with_colors=0;
    if((show_gasstarxray=='star') or (show_gasstarxray=='stars') or (show_gasstarxray=='st')): do_stars=1;
    
    
    ## check whether the data needs to be pulled up, or if its being passed by the calling routine
    if(input_data_is_sent_directly==0):
        ## read in snapshot data and do centering 
        ptypes=[2,3,4]; ## stars in non-cosmological snapshot
        if (cosmo==1): ptypes=[4];
        if (do_stars==0): ptypes=[0];
        m_all, x_all, y_all, z_all, c_all, h_all, zm_all, time = load_snapshot_brick(ptypes, \
            snapdir, snapnum, h0=use_h0, cosmological=cosmo, filename_prefix=sdir, \
            skip_bh=cosmo, do_xray=do_xray, four_char=0, use_rundir=0, full_gas_set=0);
        h_all *= 1.25;

        if ((do_stars==1) & (threecolor==1)): ## will need gas info to process attenuation
            gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_sfr, gas_lxray, gas_metallicity, gas_mass, gas_x, gas_y, gas_z, time = \
             load_snapshot_brick([0], snapdir, snapnum, h0=use_h0, cosmological=cosmo, full_gas_set=1);
            ## don't allow hot gas to have dust
            gas_temp = gadget.gas_temperature(gas_u,gas_nume); gas_metallicity[gas_temp > 1.0e6] = 0.0; 
    
        if (center[0] == 0.):
            if (cosmo==1): 
                center=gadget.calculate_zoom_center(snapdir,snapnum);
            else:
                if (center_on_bh==1):
                    pbh=gadget.readsnap(snapdir,snapnum,5,h0=use_h0,cosmological=0,skip_bh=0);
                    pos=pbh['p']; center=[pos[0,0],pos[0,1],pos[0,2]];
                if (center_on_com==1):
                    center=[np.median(x_all),np.median(y_all),np.median(z_all)];
        center=np.array(center)+np.array(offset)
        x_all-=center[0]; y_all-=center[1]; z_all-=center[2];
        print('center at ',center)


    ## rotate and re-project coordinate frame 
    if (center_is_camera_position==0):
        x=1.*x_all; y=1.*y_all; z=1.*z_all; 
        #x,y,z=coordinates_rotate(x_all,y_all,z_all,theta,phi,coordinates_cylindrical=coordinates_cylindrical);
        if ((do_stars==1) & (threecolor==1)):
            gas_x-=center[0]; gas_y-=center[1]; gas_z-=center[2];
        gx=gas_x; gy=gas_y; gz=gas_z;
        #gx,gy,gz=coordinates_rotate(gas_x,gas_y,gas_z,theta,phi,coordinates_cylindrical=coordinates_cylindrical);
    else:
        x=x_all; y=y_all; z=z_all; 
        if ((do_stars==1) & (threecolor==1)): gx=gas_x; gy=gas_y; gz=gas_z;
    
    x_tmp=1.*x; y_tmp=1.*y; z_tmp=1.*z;
    x=projection_vector_x[0]*x_tmp + projection_vector_x[1]*y_tmp + projection_vector_x[2]*z_tmp
    y=projection_vector_y[0]*x_tmp + projection_vector_y[1]*y_tmp + projection_vector_y[2]*z_tmp
    z=projection_vector_z[0]*x_tmp + projection_vector_z[1]*y_tmp + projection_vector_z[2]*z_tmp
    if ((do_stars==1) & (threecolor==1)): 
        x_tmp=1.*gx; y_tmp=1.*gy; z_tmp=1.*gz;
        gx=projection_vector_x[0]*x_tmp + projection_vector_x[1]*y_tmp + projection_vector_x[2]*z_tmp
        gy=projection_vector_y[0]*x_tmp + projection_vector_y[1]*y_tmp + projection_vector_y[2]*z_tmp
        gz=projection_vector_z[0]*x_tmp + projection_vector_z[1]*y_tmp + projection_vector_z[2]*z_tmp

    if(center_is_camera_position==0):
        x,y,z=coordinates_rotate(x,y,z,theta,phi,coordinates_cylindrical=coordinates_cylindrical);
        if ((do_stars==1) & (threecolor==1)):
            gx,gy,gz=coordinates_rotate(gx,gy,gz,theta,phi,coordinates_cylindrical=coordinates_cylindrical);



    ## set dynamic ranges of image
    temp_max=1.0e6; temp_min=1.0e3; ## max/min gas temperature
    if (do_stars==1): temp_max=50.; temp_min=0.01; ## max/min stellar age
    if (set_temp_max != 0): temp_max=set_temp_max
    if (set_temp_min != 0): temp_min=set_temp_min	
    xr=xrange; yr=xr; zr=0;
    if (checklen(yrange)>1): yr=yrange;
    if (checklen(zrange)>1): zr=zrange;
    scale=0.5*np.max([xr[1]-xr[0],yr[1]-yr[0]]); zr=np.array([-1.,1.])*scale;
    if(maxden==0.):
        maxden = 0.1 * 1.e11/1.0e10 * (2.0/scale)**(0.3) * maxden_rescale;
    if ((threecolor==1) & (do_stars)):
        if(dynrange==1.0e5): # reset if it's the default value 
            dynrange=1.0e2;
            if (nasa_colors==1): dynrange=1.0e4;
            if (nasa_colors==1): maxden *= 30.;


    ## now check if we're projecting to a camera (instead of a fixed-plane)
    xlen=0.5*(xr[1]-xr[0]); ylen=0.5*(yr[1]-yr[0]); 
    xr_0=xr; yr_0=yr;
    if (project_to_camera==1):  
        ## determine the angular opening and camera distance given physical x/y range:
        ## note: camera_opening_angle is the half-angle
        xr=np.array([-1.,1.])*np.tan(camera_opening_angle*math.pi/180.); yr=xr*ylen/xlen; ## n- degree opening angle
        ## use this to position camera
        c_dist = xlen/xr[1]; 
        camera_direction /= np.sqrt(camera_direction[0]**2.+camera_direction[1]**2.+camera_direction[2]**2.);
        camera_pos = -c_dist * camera_direction;
        ## and determine z-depth of image
        cpmax=np.sqrt(camera_pos[0]**2.+camera_pos[1]**2.+camera_pos[2]**2.); zrm=np.sqrt(zr[0]**2.+zr[1]**2.);
        ## clip z_range for speed : 
        zrm=5.0*zrm; 
        #zrm=2.5*zrm
        zrm=np.sqrt(zrm*zrm+cpmax*cpmax); zr=[-zrm,0.]
        ## correct if the center given is the camera position:
        if (center_is_camera_position==1): camera_pos=[0.,0.,0.];

        ## now do the actual projection into camera coordinates:        
        #r_pre_camera = np.sqrt(x*x+y*y+z*z)
        x,y,z,rc=coordinates_project_to_camera(x,y,z,camera_pos=camera_pos,camera_dir=camera_direction);
        h_all*=rc; ## correct for size re-scaling

        #TEMPORARY_FOR_SPIN_MOVIE_ONLY = spin_movie
        #if(TEMPORARY_FOR_SPIN_MOVIE_ONLY==1):
        #    h_max_camera = 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.) # too 'grey'? (copy 3, nothing else)
        #    h_max_camera = np.maximum(0.005+0.*z , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.)) # copy 2 
        #    h_max_camera = np.maximum(0.00707+0.*z , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.)) # copy 6
        #    h_max_camera = 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.) # good foreground; copy 7
            #h_max_camera = np.minimum(0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.), 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**0.5)) # copy 4
            #h_max_camera *= 0.7 * np.minimum((1.+(c_all/2.)**2.),(1.+(c_all/5.0)**0.75)) # copy 5 (with copy 4 line also)
            #h_max_camera *= 1. + (r_pre_camera/10.)**0.5
            #h_max_camera = np.maximum(0.005*(1.+3.0/np.abs(z/(0.1*c_dist))) , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.)) # copy 8 (1of2)
            #h_max_camera = np.maximum(0.0071*(1.+1.5/np.abs(z/(0.1*c_dist))**2.) , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**5.)) # copy 9 (1of2)
            #h_max_camera = np.maximum(0.0071*(1.+1.5/np.abs(z/(0.1*c_dist))**2.) , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**2.)) # copy 9 good except too 'sharp' a transition to smooth for later snaps (b/c of ~5 power); lower here (copy 10)
            #h_max_camera = np.maximum(0.0071/((1.+scale/2.)**0.125)*(1.+1.5/np.abs(z/(0.1*c_dist))**2.) , 0.1 * 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**2.)) # copy 9 good except too 'sharp' a transition to smooth for later snaps (b/c of ~5 power); lower here (copy 11)
            #h_max_camera = (0.0071/((1.+scale/2.)**0.125)*(0.8+0.7/np.abs(z/(0.1*c_dist))**2.)) * (1. + (np.abs(z/(0.5*c_dist)))*2.5) # good foreground; copy 7
        #    h_max_camera = 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**2.) # good foreground; copy 7_pt5 set
        #    h_all = np.minimum(h_all,h_max_camera)
        #    m_all*=1./(z*z + h_all*h_all + (0.*0.25*c_dist)**2.); ## makes into 'fluxes' (dimmer when further away)
            #m_all *= (1. + ((0.02/(0.0001+h_all))**2.) / (1. + np.abs(z/(0.1*c_dist))**2. ));  # copy 8 (2of2)
            #m_all *= (1. + ((0.02/(0.0001+h_all))**2.) / (1. + np.abs(z/(0.1*c_dist))**2. ));  # copy 9/10/11 (2of2)
            #m_all *= (1. + ((0.02/(0.0001+h_all))**2.) / (1. + np.abs(z/(0.1*c_dist))**2. ));  # copy 9/10/11 (2of2)
        #else:

        #h_max_camera = 0.02 * (h_all/0.02)**(0.5)
        #h_max_camera = 0.01 * (1. + (np.abs(z/(0.5*c_dist)))**(-2.)) # good foreground; copy 7_pt5 set
        #h_all = h_all * (h_all / 0.01)**(0.5)
        h_all_lim = 0.01
        biggie=np.where(h_all > h_all_lim); h_all[biggie] *= (h_all[biggie]/h_all_lim)**(0.0)
        #h_all = h_all * (h_all / 0.01)**(-0.5)
        h_max_camera = 1.*h_all
        #biggie = np.where(h_all > h_max_camera)
        #m_all[biggie] *= (h_max_camera[biggie]/h_all[biggie])**(0.5)
        h_all = np.minimum(h_all,h_max_camera)
        

        m_all*=1./(z*z + h_all*h_all + (0.25*c_dist)**2.); ## makes into 'fluxes' (dimmer when further away)
        ## clip particles too close to camera to prevent artifacts
        z_clip = -c_dist/10.0
        m_all[z >= z_clip] = 0.0
        ## also need to correct gas properties if using it for attenuation
        if ((threecolor==1) & (do_stars==1)):
            gx,gy,gz,grc=coordinates_project_to_camera(gx,gy,gz,camera_pos=camera_pos,camera_dir=camera_direction);
            gas_hsml*=grc; gas_mass*=1./(gz*gz + gas_hsml*gas_hsml + (0.25*c_dist)**2.); ## this conserves surface density, so attenuation units are fine
            gas_mass[gz >= z_clip] = 0.0 ## clipping to prevent artifacts
            biggie=np.where(gas_hsml > h_all_lim); gas_hsml[biggie] *= (gas_hsml[biggie]/h_all_lim)**(0.0)
        
        #rotate the image; assumes angle is specified in degrees
        angle_to_rotate_image*=math.pi/180.
        xnew = x*np.cos(angle_to_rotate_image) - y*np.sin(angle_to_rotate_image)
        ynew = x*np.sin(angle_to_rotate_image) + y*np.cos(angle_to_rotate_image)
        x = xnew
        y = ynew
        if ((threecolor==1) & (do_stars==1)):
            gxnew = gx*np.cos(angle_to_rotate_image) - gy*np.sin(angle_to_rotate_image)
            gynew = gx*np.sin(angle_to_rotate_image) + gy*np.cos(angle_to_rotate_image)
            gx = gxnew
            gy = gynew

    #plotting_setup_junk
    plt.close('all')
    format = '.pdf' #format = '.png' ## '.ps','.pdf','.png' also work well
    axis_ratio = ylen/xlen
    #fig=plt.figure(frameon=False,figsize=(10.,10.*axis_ratio),dpi=100)
    ## for Chris's high-res images:
    fig=plt.figure(frameon=False,dpi=pixels,figsize=(1.,1.*axis_ratio))
    fig.set_size_inches(set_size_inches,set_size_inches*axis_ratio)
    ax_fig=plt.Axes(fig,[0.,0.,1.,1.*axis_ratio],clip_on=True)
    if(axis_enabled==False): ax_fig.set_axis_off()
    frame1=plt.gca()
    if(axis_enabled==False): frame1.set_axis_off()
    plt.xlabel(''); plt.ylabel('')
    frame1.axes.xaxis.set_ticklabels([]) #no tick names
    frame1.axes.yaxis.set_ticklabels([]) #no tick names
    frame1.axes.get_xaxis().set_ticks([]) #no ticks
    frame1.axes.get_yaxis().set_ticks([]) #no ticks
    ax_fig.axes.xaxis.set_ticklabels([]) #no tick names
    ax_fig.axes.yaxis.set_ticklabels([]) #no tick names
    ax_fig.axes.get_xaxis().set_ticks([]) #no ticks
    ax_fig.axes.get_yaxis().set_ticks([]) #no ticks
    if(axis_enabled==True): fig.add_axes(ax_fig)
    if(axis_enabled==False): ax_fig=frame1

	  
	## trim particles to those inside/near the plotting region
    sidebuf=10.; 
    if (cosmo==1): sidebuf=1.15
    x=np.array(x,dtype='f'); y=np.array(y,dtype='f'); z=np.array(z,dtype='f'); 
    m_all=np.array(m_all,dtype='f'); h_all=np.array(h_all,dtype='f'); 
    c_all=np.array(c_all,dtype='f'); 
    ok =    ok_scan(x,xmax=np.max(np.fabs(xr))*sidebuf) & ok_scan(y,xmax=np.max(np.fabs(yr))*sidebuf) & \
            ok_scan(z,xmax=np.max(np.fabs(zr))*sidebuf) & ok_scan(m_all,pos=1,xmax=1.0e40) & \
            ok_scan(h_all,pos=1) & ok_scan(c_all,pos=1,xmax=1.0e40)
    weights = m_all; color_weights = c_all; ## recall, gadget masses in units of 1.0d10
    print("Mean h_all = ",h_all.mean())
    print("Median h_all = ",np.median(h_all))


    TEMPORARY_FOR_THORSTEN_REVIEW_ONLY = 0
    if(TEMPORARY_FOR_THORSTEN_REVIEW_ONLY==1):
        #xxx=1.*x; x=1.*y; y=1.*xxx;
        if ((threecolor==1) & (do_stars==1)):
            #gxxx=1.*gx; gx=1.*gy; gy=1.*gxxx;
            h_all *= 0.8
            #c_all *= (c_all / 1.0);
            gas_metallicity *= 1.e-5;
            c_all *= (c_all / 1.0)**0.5; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.25; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.5; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.25; 
            use_old_extinction_routine = 0;

    TEMPORARY_FOR_SPIN_MOVIE_ONLY = spin_movie
    if(TEMPORARY_FOR_SPIN_MOVIE_ONLY==1):
        if ((threecolor==1) & (do_stars==1)):
            h_all *= ((scale/30.)**0.25) * (1.0 + (h_all / 0.1)**(0.5)) / (1. + (h_all / 1.0)**(0.5))
            #c_all *= (c_all / 1.0)**0.25; 
            c_all *= (c_all / 1.0)**0.125; 
            #h_all *= h_rescale_factor
            ## limit softenings to a maximum fraction of FOV to avoid foreground fuzziness/blobs
            #h_all = np.minimum(h_all,h_max) # CCH
        else:
            h_all *= (((30.*0.+1.*scale)/30.)**0.25) * (0.8 + (h_all / 0.1)**(0.5)) / (1. + (h_all / 1.0)**(0.5))

    ## adjust softenings if appropriate keywords are set
    h_all *= h_rescale_factor
    if (h_max > 0):
        h_all = np.minimum(h_all,h_max)

    ## alright, now ready for the main plot construction:
    if (threecolor==1):
        if (do_stars==1):
            ## making a mock 3-color composite image here ::
            
            ## first grab the appropriate, attenuated stellar luminosities
            ## for now not including BHs, so initialize some appropriate dummies:
            bh_pos=[0,0,0]; bh_luminosity=0.; include_bh=0;
            ## and some limits for the sake of integration:
            SizeMax=20.*np.max([np.fabs(xr[1]-xr[0]),np.fabs(yr[1]-yr[0])]);
            SizeMin=np.min([np.fabs(xr[1]-xr[0]),np.fabs(yr[1]-yr[0])])/250.;
 
            ## now set up the main variables:
            star_pos=np.zeros((3,checklen(x[ok]))); 
            star_pos[0,:]=x[ok]; star_pos[1,:]=y[ok]; star_pos[2,:]=z[ok];
            stellar_age=np.maximum(c_all,min_stellar_age);
            stellar_metallicity=zm_all; stellar_mass=m_all;
            print('Min stellar age = ',stellar_age.min())
            print('Mean stellar age = ',stellar_age.mean())
            print('Mean stellar metallicity = ',stellar_metallicity.mean())
            ## gas will get clipped in attenuation routine, don't need to worry about it here
            gas_pos=np.zeros((3,checklen(gx))); 
            gas_pos[0,:]=gx; gas_pos[1,:]=gy; gas_pos[2,:]=gz;

            if ((add_gas_layer_to_image=='')==False):
                gas_wt = 0.*gx
                maxden_for_layer=1.*maxden; dynrange_for_layer=1.*dynrange
                set_percent_maxden_layer=0.; set_percent_minden_layer=0.;
                gas_hsml_for_extra_layer=gas_hsml
                if (add_gas_layer_to_image=='Halpha'): 
                    gas_wt = gas_mass * gas_rho*gas_nume*(1.-gas_numh)*(gadget.gas_temperature(gas_u,gas_nume)**(-0.75))/(gadget.gas_mu(gas_nume)**2.)
                    maxden_for_layer *= 3.e4;
                    #dynrange_for_layer *= 1.5;
                    dynrange_for_layer *= 0.7;
                    dynrange_for_layer *= 2.0;
                    #gas_hsml_for_extra_layer *= 1.1;
                if (add_gas_layer_to_image=='CO'): 
                    gas_tmp = gadget.gas_temperature(gas_u,gas_nume)
                    n_tmp = gas_rho * 176.2; # assuming mean molec weight of 2.3 for dense gas
                    gas_wt = gas_mass * np.exp(-(gas_tmp/8000. + 10./n_tmp));
                    maxden_for_layer *= 3.e4;
                    dynrange_for_layer *= 0.7e4;
                if (add_gas_layer_to_image=='SFR'): 
                    gas_wt = gas_sfr
                    set_percent_maxden_layer=0.9999;
                    set_percent_minden_layer=0.01;
                if (add_gas_layer_to_image=='Xray'): 
                    gas_wt = gas_lxray
                    maxden_for_layer *= 0.6e3;
                    dynrange_for_layer *= 0.1;
                    gas_hsml_for_extra_layer *= 1.2;
                if (add_gas_layer_to_image=='Zmetal'): 
                    gas_wt = gas_mass*gas_metallicity
                    set_percent_maxden_layer=0.9999;
                    set_percent_minden_layer=0.01;
                gas_wt /= np.sum(gas_wt)
                                    
                massmap_gas_extra_layer,image_singledepth_extra_layer = \
                cmakepic.simple_makepic(gx,gy,weights=gas_wt,hsml=gas_hsml_for_extra_layer,\
                    xrange=xr,yrange=yr,
                    set_dynrng=dynrange_for_layer,set_maxden=maxden_for_layer,
                    set_percent_maxden=set_percent_maxden_layer,set_percent_minden=set_percent_minden_layer, 
                    color_temperature=0,pixels=pixels,invert_colorscale=1-invert_colors);


            if (use_old_extinction_routine == 0):
                ##
                ## this is the newer 'single call' attenuation and ray-tracing package: 
                ##  slightly more accurate in diffuse regions, more stable behavior
                ##

                if (x[ok].size <= 3):
                    ## we got a bad return from the previous routine, initialize blank arrays
                    out_gas=out_u=out_g=out_r=np.zeros((pixels,pixels))
                    image24=massmap=np.zeros((pixels,pixels,3))
                else:
                    ## actually call the ray-trace:
                    out_gas,out_u,out_g,out_r = rayproj.stellar_raytrace( BAND_IDS, \
                        x[ok], y[ok], z[ok], \
                        stellar_mass[ok], stellar_age[ok], stellar_metallicity[ok], h_all[ok], \
                        gx, gy, gz, gas_mass, \
                        gas_metallicity*dust_to_gas_ratio_rescale, gas_hsml, \
                        xrange=xr, yrange=yr, zrange=zr, pixels=pixels, \
                        #ADD_BASE_METALLICITY=1.0*0.02, ADD_BASE_AGE=0.003, 
                        ADD_BASE_METALLICITY=1.0e-6*0.02, ADD_BASE_AGE=0.0003, 
                        #ADD_BASE_METALLICITY=0.01*0.02, ADD_BASE_AGE=0.0003, 
                        IMF_SALPETER=0, IMF_CHABRIER=1 );
                
                    if(np.array(out_gas).size<=1):
                        ## we got a bad return from the previous routine, initialize blank arrays
                        out_gas=out_u=out_g=out_r=np.zeros((pixels,pixels))
                        image24=massmap=np.zeros((pixels,pixels,3))
                    else:
                        ## make the resulting maps into an image
                        image24, massmap, image_max, image_min = \
                            makethreepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, \
                            maxden=maxden,dynrange=dynrange,pixels=pixels, \
                            color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors, \
                            set_f_saturated=set_f_saturated);  
                
            else: ## use_old_extinction_routine==1; can be useful for certain types of images
                ## (for example, this can allow for a scattered fraction, the other does not currently)

                ## call routine to get the post-attenuation band-specific luminosities
                lum_noatten, lum_losNH = \
                getlum.get_attenuated_stellar_luminosities( BAND_IDS, star_pos, gas_pos, bh_pos, \
                    stellar_age[ok], stellar_metallicity[ok], stellar_mass[ok], \
                    gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, \
                    gas_metallicity*dust_to_gas_ratio_rescale, gas_mass, \
                    bh_luminosity, \
                    xrange=xr, yrange=yr, zrange=zr, \
                    INCLUDE_BH=include_bh, SKIP_ATTENUATION=0, 
                    #ADD_BASE_METALLICITY=1.0, ADD_BASE_AGE=0., 
                    ADD_BASE_METALLICITY=1.0e-2, ADD_BASE_AGE=0., 
                    IMF_SALPETER=0, IMF_CHABRIER=1, \
                    MIN_CELL_SIZE=SizeMin, OUTER_RANGE_OF_INT=SizeMax, \
                    SCATTERED_FRACTION=scattered_fraction, \
                    REDDENING_SMC=1, REDDENING_LMC=0, REDDENING_MW=0, \
                    AGN_MARCONI=0, AGN_HRH=1, AGN_RICHARDS=0, AGN_SDSS=0 );
                
                ## now pass these into the 3-band image routine
                image24, massmap = \
                makethreepic.make_threeband_image(x[ok],y[ok],lum_losNH,hsml=h_all[ok],\
                    xrange=xr,yrange=yr,maxden=maxden,dynrange=dynrange,pixels=pixels, \
                    color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors );           

        else: 
            ## 
            ## threecolor==1, but do_stars==0, so doing a multi-pass gas image
            ## 
            ##  -- many experiments here: doing gas isosurfaces with broad kernels
            ##       and overlaying a custom set of color tables after the fact seems
            ##       best. adjust opacity (kappa_units) and kernel width as needed. 
            ##       also, can set 'dynrange=0' to automatically produce a given 
            ##       fraction of saturated/unsaturated pixels
            ##
            gas_map_temperature_cuts=np.array([300., 2.0e4, 3.0e5 ])
            kernel_widths=np.array([0.8,0.3,0.6])
            #kernel_widths=np.array([0.7,0.3,0.7])

            # layer below is for zoom_dw large-scale GAS sims ONLY
            #gas_map_temperature_cuts=np.array([300., 1.0e4, 1.0e5 ])
            #kernel_widths=np.array([0.5,0.25,0.6])

            out_gas,out_u,out_g,out_r = rayproj.gas_raytrace_temperature( \
                gas_map_temperature_cuts, \
                x[ok], y[ok], z[ok], color_weights[ok], weights[ok], h_all[ok], \
                xrange=xr, yrange=yr, zrange=zr, pixels=pixels, \
                isosurfaces = 1, kernel_width=kernel_widths, \
                add_temperature_weights = 0 , KAPPA_UNITS = 2.0885*np.array([1.1,2.0,1.5]));
                #add_temperature_weights = 0 , KAPPA_UNITS = 2.0885*np.array([4.1,2.0,2.0]));
                
            image24, massmap, image_max, image_min = \
                makethreepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, \
                maxden=maxden,dynrange=dynrange,pixels=pixels, \
                color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors);

            image24 = makethreepic.layer_band_images(image24, massmap);                

    else:   ## threecolor==0    ( not a three-color image )
    
        if (do_with_colors==1):
            ##
            ## here we have a simple two-pass image where surface brightness is 
            ##   luminosity and a 'temperature' is color-encoded
            ##
            if (log_temp_wt==1): 
                color_weights=np.log10(color_weights); temp_min=np.log10(temp_min); temp_max=np.log10(temp_max);

            massmap_singlelayer,pic_singlelayer,massmap,image24 = \
            cmakepic.simple_makepic(x[ok],y[ok],weights=weights[ok],hsml=h_all[ok],\
                xrange=xr,yrange=yr,
                set_dynrng=dynrange,set_maxden=maxden,
                set_percent_maxden=set_percent_maxden,set_percent_minden=set_percent_minden, 
                color_temperature=1,temp_weights=color_weights[ok],
                pixels=pixels,invert_colorscale=invert_colors,
                set_temp_max=temp_max,set_temp_min=temp_min);

        else: ## single color-scale image
            ##
            ## and here, a single-color scale image for projected density of the quantity
            ##
            massmap,image_singledepth = \
            cmakepic.simple_makepic(x[ok],y[ok],weights=weights[ok],hsml=h_all[ok],\
                xrange=xr,yrange=yr,
                set_dynrng=dynrange,set_maxden=maxden,
                set_percent_maxden=set_percent_maxden,set_percent_minden=set_percent_minden, 
                color_temperature=0,pixels=pixels,invert_colorscale=1-invert_colors);

            ## convert to actual image using a colortable:
            my_cmap=matplotlib.cm.get_cmap('hot');
            image24 = my_cmap(image_singledepth/255.);

    ##
    ## whichever sub-routine you went to, you should now have a massmap (or set of them) 
    ##   and processed rgb 'image24' (which should always be re-makable from the massmaps)
    ##
    ## optionally can have great fun with some lighting effects to give depth:
    ##   (looks great, with some tuning)
    ##
    if (include_lighting==1):
        #light = matplotlib.colors.LightSource(azdeg=0,altdeg=65)
        light = viscolors.CustomLightSource(azdeg=0,altdeg=65)
        if (len(massmap.shape)>2):
            ## do some clipping to regulate the lighting:
            elevation = massmap.sum(axis=2)
            minden = maxden / dynrange
            elevation = (elevation - minden) / (maxden - minden)
            elevation[elevation < 0.] = 0.
            elevation[elevation > 1.] = 1.
            elevation *= maxden
            grad_max = maxden / 5.
            grad_max = maxden / 6.
            #image24_lit = light.shade_rgb(image24, massmap.sum(axis=2))
            image24_lit = light.shade_rgb(image24, elevation, vmin=-grad_max, vmax=grad_max)
        else:
            image24_lit = light.shade(image24, massmap)
        image24 = image24_lit
    
    plt.imshow(image24,origin='lower',interpolation='bicubic',aspect='normal');

    if ((add_gas_layer_to_image=='')==False):
        viscolors.load_my_custom_color_tables();
        plt.imshow(image_singledepth_extra_layer,origin='lower',interpolation='bicubic',aspect='normal',
            cmap=set_added_layer_ctable,alpha=set_added_layer_alpha)

    ## slap on the figure labels/scale bars
    if(show_scale_label==1): overlay_scale_label(xr_0,yr_0,ax_fig,c0='w')
    if(show_time_label==1): overlay_time_label(time,\
        ax_fig,cosmo=cosmo,c0='w',n_sig_add=0,tunit_suffix='Gyr')    
    
    filename=fname_base
    if(filename_set_manually!=''): filename=filename_set_manually
    plt.subplots_adjust(left=0.0,bottom=0.0,right=1-0.0,top=1-0.0,hspace=0.0,wspace=0.0);
    frame1.axes.xaxis.set_ticklabels([]) #no tick names
    frame1.axes.yaxis.set_ticklabels([]) #no tick names
    frame1.axes.get_xaxis().set_ticks([]) #no ticks
    frame1.axes.get_yaxis().set_ticks([]) #no ticks
    plt.tight_layout(pad=0,w_pad=0,h_pad=0)
    plt.savefig(filename+format,dpi=pixels,bbox_inches='tight',pad_inches=0)
    plt.close('all')

    ## save them:
    outfi = h5py.File(filename+'.dat','w')
    dset_im = outfi.create_dataset('image24',data=image24)
    dset_mm = outfi.create_dataset('massmap',data=massmap)
    dset_xr = outfi.create_dataset('xrange',data=xr_0)
    dset_yr = outfi.create_dataset('yrange',data=yr_0)
    dset_time = outfi.create_dataset('time',data=time)
    dset_cosmo = outfi.create_dataset('cosmological',data=cosmo)
    if ((add_gas_layer_to_image=='')==False):
        dset_im2 = outfi.create_dataset('image24_extralayer',data=image_singledepth_extra_layer)
    else:
        dset_im2 = outfi.create_dataset('image24_extralayer',data=np.array([-1]))
    outfi.close()
    ## (to read) :::
    #infiname = fname_base+'.dat'
    #infi=h5py.File(infiname,'r')
    #image24 = np.array(infi["image24"])
    #massmap = np.array(infi["massmap"])
    #infi.close()
    
    ## and return the arrays to the parent routine!
    return image24, massmap;


def xcompile(x_grid_previous,x_final,dlnx_previous,n):
    x_0 = x_grid_previous[-1]
    x_f = x_final
    dx_desired = np.log(x_f/x_0) / n
    lnx_grid_0 = np.log(x_grid_previous)
    lnx_new, dlnx_new = tcompile(lnx_grid_0,dx_desired,dlnx_previous,n)
    return np.exp(lnx_new), dlnx_new
    #return np.concatenate([x00_grid, x00_i + (x00_f-x00_i)/n*np.arange(0,n)])
    #return np.concatenate([x00_grid, x00_i * np.exp( np.log(x00_f/x00_i)/n*np.arange(0,n) )])

def tcompile(t_grid,dt_desired,dt_previous,n):
    f_transition = 0.1
    dt_0 = dt_previous
    dt_1 = (dt_desired - 0.5*f_transition*dt_previous)/(1.0 - 0.5*f_transition)
    dt = np.zeros(n)
    ng = 1.*np.arange(1,n+1)
    ok = (ng <= f_transition*n)
    dt[ok] = dt_0 + (dt_1-dt_0) * (1.*ng[ok])/(f_transition*n)
    ok = (ng > f_transition*n)
    dt[ok] = dt_1
    tg = np.cumsum(dt)
    t_grid_new = np.concatenate([t_grid, t_grid[-1] + tg])
    return t_grid_new, dt_1



def frame_ext(snum,four_char=1):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;

