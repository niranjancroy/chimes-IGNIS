import numpy as np
import h5py as h5py
import os.path
import math
import matplotlib
#matplotlib.use('Agg') ## this calls matplotlib without a X-windows GUI
import matplotlib.pylab as pylab
import matplotlib.pyplot as plot
import scipy.misc
import scipy.interpolate as interpolate
import scipy.optimize as optimize
#import astropy.cosmology as cosmo
from astropy.cosmology import WMAP9 as cosmo
import fit_functions as fitfun
#cosmo.set_current('WMAP9')
import gadget
import colors_sps.colors_table as sps
import halo_finder
#import atpy
import cosmology as mycosmo
import pfh_utils as util
import scipy.special as scifn
#import asciitable
from astropy.io import ascii
import astropy
import scipy.stats


def master():
	sdir='.'; 
	snum=440;
	cen=[3306.56,3029.35,3235.48];
	
	sdir_m='../../zooms/'
	sdir_m='../zooms/'
	
	sdir='m12_mr'
	snum=364
	snum=441
	#snum=429
	#snum=418 # good 'before merger' snapshot
	#sdir='hires_bh'
	#snum=441
	#sdir='m13_ics'
	#snum=56
	sdir='m11_hr'
	snum=316 #440
	snum=440
	sdir='m11_lr'
	snum=440
	#sdir='../zooms/m10_dw'
	#snum=440
	
	sdir=sdir_m+sdir

	
	cen=calculate_zoom_center(sdir,snum);
	print('found center at ',cen)
	
	#quicklook(sdir,snum,cen,5.);
	#dm_profile(sdir,snum,cen);
	#vc_profile(sdir,snum,cen);
	#sfr_time(sdir,snum,cen);
	#sb_profile(sdir,snum,cen);
	ang_mom(sdir,snum,cen);



def local_get_halo_id_filename(snapdir, snum, full_halo_results=0, four_char=0):
    ## (first parse the directory names (for file-naming conventions))
    s0=snapdir.split("/"); snapdir_specific=s0[len(s0)-1]; n_s0=1;
    if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2]; n_s0=2;
    ss=halo_finder.snap_ext(snum,four_char=four_char)
    fname_base='./halos/'+snapdir_specific+ss
    
    ext = '.dmhalo_ids'
    if (full_halo_results==1): ext = '.dmhalo_prop'
    
    return fname_base+ext

"""
 Mgal-Mhalo plots:
    - relation at z=0 with 'all' galaxies
        - panel with 'efficiency' also? (so can see how low it goes?)
    - zoom in on that plot: compare w. other sims (low-mass and high-mass ranges)
    - numerical effects: show results w. different res, numerical parameters, etc.
    - physical effects: show results w. different feedback models
        - related: show sfr vs time? worth doing, even if will be repeated in SFH ppr? 
    - evolution w. redshift (panels at different z)
    - scatter 
        - show distrib?: is it log-normal? 
        - scatter (in dex) in bins vs. gal or halo mass 
        - possibly vs. z too - dift colors on same plot, or plot with mean trend?
    - for dwarfs: implications vs 'too big to fail'? (own paper?)
        - (would be mg-mh plot from MBK+Bullock to compare)
    - MF: convolve our pred at each z, with scatter, with halo MF to get galaxy MF
        - (truncated at most massive gal we simulate)
        - (compare to obs at each z; redundant info, but worth it)
    - could also predict clustering for each Mgal/z? Or just fitted Mgal-Mhalo to quote?
    - want to include K-S here? also 'self-regulation' of star formation....?
        (but it does get a bit into the SFH aspect of things...)
    - physics? outflow and inflow rates???
        
        
 SFH plots:
    - example SFR vs. time (m12 sim): possibly series of panels with different masses?
        - compare no fb, w. fb, w. sub-grid fb runs (one figure)
        - 'typical' SFH in another (w panel vs mass, so can show some features)
    - SFR 'main sequence' 
        - z=0, and same fig or others with panels at different z
    - numerical parameters: sfr sequence vs. knobs (or specific examples?)
    - physical parameters: sfr sequence vs. dif't feedback
    - convolve w. halo MF -> madau plot vs. redshift (do & don't cut off at max masses we model)
    - 'burstiness': (show?) SFR with dif't time-averaging: 
        - scatter in SFR vs. time vs. dift time-avg window size
           (do for several sims, plot vs. dt if all same, unless systematic dependence on 
            mass/z, then plot in panels each with that, for comparison)
    - 'tau model' properties:
        - 'effective age' vs. gal mass (w and w/o free metallicity)
        - 'decline time' vs. mass (show for dif't fits; compare to noeske, etc.)
            - show vs. z: at high-z, show preference for rising SFHs (negative tau)
    - Kennicutt relations:
        - all gals at z=0, panels vs. redshift
        - at z~0, do 'resolved' in regions inside of galaxies?
    - compare to halo growth rates * universal baryon fraction:
        - show SFR vs time, vs this, and 'corrected' version (forcing to live on mean mg-mh relation)
            (use actual halo growth rate in sim, so can compare acc. rate)
        - show same for 'main sequence' (panels)
        - show Mdot 'into galaxy'? vs SFR or dMhalo/dt? (this gets a bit into dusan's territory, 
            don't want to go too far into what determines galaxy accretion rates....)
        - compare inflow and outflow rates? 'whats needed' to get low SFRs?
            (also steps a bit on other projects)
    - what causes burstiness? 
        - 'regular' burstiness (scatter within main sequence): 
            - fb effects? (compare no-fb or decoupled winds model)
            - disk instabilities (measure mode amplitudes in disk at 'peak' vs. 'trough' times, 
                compare to see if this correlates with bursts)
            - mergers? (for strongest bursts?):
                - definitely emphasize see stronger FX here than in 'sub-grid' FB models
                - correlate peaks w. companions? correlate companions w. peaks? 
                    ( show SFH when merging & non-merging times marked ? )

"""
def plot_mgal_mhalo():
    import halo_finder as hfind
    snapdir_m='/Users/phopkins/Documents/work/plots/zooms/'
    pylab.close('all')
    
    SHOW_MGAL_MHALO_Z0 = 0
    SHOW_MGAL_MHALO_ALLZ = 1
    SHOW_MGAL_MHALO_NUMERICS = 0
    SHOW_MGAL_MHALO_ALTSIMS = 0
    SHOW_MGAL_SCATTER = 0
    SHOW_SFR_MAINSEQUENCE = 0
    SHOW_SSFR = 0
    SHOW_MGAL_METALLICITY = 0
    SHOW_KSLAW = 0

    z_vec = [0.0]
    if (SHOW_MGAL_MHALO_Z0): z_vec=[0.0]
    if (SHOW_MGAL_MHALO_NUMERICS): z_vec=[0.0, 0.05]
    if (SHOW_MGAL_MHALO_ALLZ): z_vec = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0]
    if (SHOW_MGAL_SCATTER): z_vec = [0.0]#, 0.5, 1.0, 2.0, 4.0, 6.0]
    if (SHOW_MGAL_MHALO_ALTSIMS): z_vec = [0.0, 2.0]
    if (SHOW_SFR_MAINSEQUENCE): z_vec = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0]
    if (SHOW_MGAL_METALLICITY): z_vec = [0.0]
    if (SHOW_SSFR): z_vec = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0]
    if (SHOW_KSLAW): z_vec = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0]

    if ((SHOW_SFR_MAINSEQUENCE) | (SHOW_SSFR)):
        plot.figure(1,figsize=(16.,6.))
        #plot.figure(1,figsize=(16.,8.))
        charsi=20.
        matplotlib.rcParams.update({'font.size':18})
    
    if (SHOW_MGAL_METALLICITY):
        plot.figure(1,figsize=(8.,6.))
        charsi=20.
        matplotlib.rcParams.update({'font.size':18})

    if (SHOW_KSLAW):
        plot.figure(1,figsize=(8.,6.))
        charsi=20.
        matplotlib.rcParams.update({'font.size':18})
    
    if(SHOW_MGAL_MHALO_ALLZ):
        #plot.figure(1,figsize=(24.,6.))
        plot.figure(1,figsize=(12.,6.))
        charsi=18.
        matplotlib.rcParams.update({'font.size':16})

    if(SHOW_MGAL_SCATTER):
        plot.figure(1,figsize=(8.,6.))
        charsi=20.
        matplotlib.rcParams.update({'font.size':18})
    
    default_B1='B1_hr_Aug27_2013_egylim'#'B1_hr_Aug8_2013_',hires_bh,B1_hr_Aug27_2013_egylim,B1_hr_Aug27_2013_egylim4
    default_m09='m09_hr_Aug8_2013_' # m09_dw,m09_hr_Aug8_2013_
    default_m10='m10_hr_Aug8_2013' #zoom_dw,'m10_Aug8_2013_'
    default_m11='m11_hr_Aug8_2013_' #m11_gas_hr,m11_gas_new
    default_m12='m12_hr_Sep10_2013_nolim'#'m12_mr_Aug8_2013',m12_mr
    default_m13='m13_mr_Aug8_2013' #m13_gas_lr


    default_B1='B1_hr_Dec5_2013_11_'#'B1_hr_Aug8_2013_',hires_bh,B1_hr_Aug27_2013_egylim,B1_hr_Aug27_2013_egylim4
    default_m09='m09_hr_Dec16_2013_' # m09_dw,m09_hr_Aug8_2013_
    default_m10lr='m10_lr_Dec9_2013_' #zoom_dw,'m10_Aug8_2013_'
    default_m10='m10_hr_Dec9_2013_' #zoom_dw,'m10_Aug8_2013_'
    default_m11lr='m11_hr_Dec9_2013_' #m11_gas_hr,m11_gas_new
    default_m11='m11_hhr_Jan9_2013_' #m11_gas_hr,m11_gas_new
    default_m12='m12qq_hr_Dec16_2013_'#'m12_mr_Aug8_2013',m12_mr
    default_m12v='m12v_mr_Dec5_2013_3_'#'m12_mr_Aug8_2013',m12_mr
    default_m13='m13_mr_Dec16_2013_' #m13_gas_lr
    default_m14='m13m14_lr_Dec9_2013_' #m13_gas_lr

    default_m10='m10v_lr_Jan2014_' #zoom_dw,'m10_Aug8_2013_'
    default_m11='m11v_lr_Jan2014_' #m11_gas_hr,m11_gas_new

    
    for i_z,z_redshift in zip(range(np.array(z_vec).size),z_vec):
        #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_hr', 'm09_dw', 'm13_gas_lr','m13_gas_lr','m14_lr','m14_lr']
        sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13]#,'m13_gas_lr','m14_lr','m14_lr']
        sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v,default_m14]
        sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v]#,default_m14]
        if (SHOW_MGAL_MHALO_Z0):
            #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_hr', 'm09_dw', 'm13_gas_lr','m13_gas_lr']#,'m14_lr','m14_lr']
            sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13]#,'m13_gas_lr']#,'m14_lr','m14_lr']
            sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v]#,default_m14]
        sdirlbl=[ 'm12v'    , 'm10' ,   'm12q'  ,   'm11',      'm09',      'm13',      'null',     'm14',  'null','null','null','null','null','null','null']
        sdirlbl=[ 'm12v'    , 'm10' ,   'm12q'  ,   'm11',      'm09',      'm13',      'm12i',     'm14',  'null','null','null','null','null','null','null']
        if (z_redshift >= 0.5):
            #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_hr', 'm09_dw', 'm13_gas_lr','m13_gas_lr','m14_lr','m14_lr']
            sdirs=[default_B1, default_m10, default_m12, default_m11, default_m09, default_m13,default_m13,'m14_lr','m14_lr']
            sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v]#,default_m14]
        if (SHOW_MGAL_MHALO_ALLZ)|(SHOW_MGAL_MHALO_ALTSIMS):
            #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_hr', 'm09_dw', 'm13_gas_lr','m13_gas_lr']#,'m14_lr','m14_lr']
            sdirs=[default_B1, default_m10, default_m12, default_m11, default_m09, default_m13]#,'m13_gas_lr']#,'m14_lr','m14_lr']
            sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v]#,default_m14]
        m_primary_list=[13.,  10.,      13.5,        12.7,            10.5,     15.,10.,                14.2,    14.9]
        m_primary_list=[13.2,  8.,      13.2,        12.0,            9.,     14.,13.3]#,                15.0,    14.9]
        contam_cuts=[ 0.1,    0.5 ,       0.5   ,   0.1     ,       0.5,     1.,1.,                 0.2, 0.2, 0.2, 0.2]
        contam_cuts=[ 0.1,    0.07 ,       0.5   ,   0.1     ,       0.5,     0.1,0.1,                 0.1]#, 0.2, 0.2, 0.2]
        if ((SHOW_SFR_MAINSEQUENCE) | (SHOW_SSFR)):
            #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_new', 'm09_dw', 'm13_gas_lr','m13_gas_lr','m14_lr','m14_lr']
            sdirs=['B1_hr_Aug8_2013_', 'zoom_dw', default_m12, default_m11, default_m09, default_m13,default_m13,'m14_lr','m14_lr']
            sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v,default_m14]
            m_primary_list=[10.,  10.,        10.,             10.,          10.0,     10.,        10.,                14.2,    14.9]
            m_primary_list=[8.2,  8.,          10.0,           10.0,            9.,     10.,       10,                13.5,    14.9]
            contam_cuts=0.+np.array([0.1,0.5 , 0.5   ,        0.1     ,       0.5,     1.,          1.,               0.01, 0.2, 0.2, 0.2])
        if (SHOW_MGAL_MHALO_Z0)|(SHOW_MGAL_MHALO_ALTSIMS):
            #m_primary_list=[12.5,  10.,      12.2,        11.,            10.8,     10.,10.,                14.2,    14.9]
            m_primary_list=[12.3,  8.,      12.7,        11.,            9.,     15.,10.,                14.2,    14.9]
            m_primary_list=[13.2,  8.,      13.5,        13.1,            9.,     15.,10.,                14.2,    14.9]
            m_primary_list=[13.2,  8.,      13.2,        11.1,            9.,     15.5,14.15,                15.7]#,    14.9]
        if (SHOW_SFR_MAINSEQUENCE==0)&(SHOW_SSFR==0)&(z_redshift >= 1.0):
            #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_hr', 'm09_dw', 'm13_gas_lr','m13_gas_mr','m14_lr','m14_lr']
            sdirs=['B1_hr_Aug8_2013_', default_m10, default_m12, default_m11, default_m09, default_m13,'m13_gas_mr','m14_lr','m14_lr']
            sdirs=['B1_hr_Aug8_2013_', default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v,default_m14]
            sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v]#,default_m14]
            if (SHOW_MGAL_MHALO_ALLZ)|(SHOW_MGAL_MHALO_ALTSIMS):
                #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_hr', 'm09_dw', 'm13_gas_lr','m13_gas_mr']#,'m14_lr','m14_lr']
                sdirs=[default_B1, default_m10, default_m12, default_m11, default_m09, default_m13]#,'m13_gas_mr']#,'m14_lr','m14_lr']
                sdirs=[default_B1, default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v]#,default_m14]#,'m13_gas_mr']#,'m14_lr','m14_lr']
            m_primary_list=[13.,  10.,      12.2,        12.,            10.,     13.5,13.5,                14.2,    14.9]
            m_primary_list=[13.,  10.,      13.5,        12.7,            10.5,     15.,10.,                14.2,    14.9]
            m_primary_list=[13.,  10.,      12.2,        12.0,            10.0,     13.9,12.5]#,                14.5,    13.]
            if (SHOW_MGAL_MHALO_ALTSIMS):
                sdirs=[default_B1, default_m10, default_m12, default_m11, default_m09, default_m13,default_m12v]#,default_m14]#,'m13_gas_mr']#,'m14_lr','m14_lr']
                m_primary_list=[13.,  10.,      12.2,        12.0,            10.0,     14.1,        12.5]#,                14.5,    13.]
            if (z_redshift >= 4.0):
                m_primary_list=[13.,  10.,      12.2,        12.0,            10.0,     13.5,12.5,                13.8,    13.]
            contam_cuts=[ 0.1,    0.5 ,       0.5   ,   0.1     ,       0.5,     1.0,0.1,                 0.2, 0.2, 0.2, 0.2]
            contam_cuts=[ 0.1,    0.07 ,       0.5   ,   0.1     ,       0.5,     0.1,0.1,                 0.2, 0.2, 0.2, 0.2]

        #snums_z0=[441,             440,      440  ,    440     ,       339,      66,66,                275,312]
        #snums_z05=[340,             340,      340  ,    340     ,    339,      57,48 ,                  275,312         ]
        #snums_z0=[440,             440,      440  ,    440     ,       440,      427,66,                275,312]
        snums_z0=[440,             440,      440  ,    440     ,       440,      427,427,                275,312]
        snums_z0=[423,             440,      440  ,    340     ,       440,      440,440,                440,440]
        snums_z0=[440,             440,      440  ,    340     ,       440,      440,440,                440,440]
        snums_z05=[340,             340,      340  ,    340     ,    340,      340,340 ,                  275,312         ]
        snums_z05=[340,             340,      340  ,    340     ,    340,      340,340 ,                  340,340         ]
        #snums_z0=[440,             440,      440  ,    440     ,       339,      427,66,                275,312]
        #snums_z05=[340,             340,      340  ,    340     ,    339,      340,48 ,                  275,312         ]
        snums_z1=[290,             290,      290  ,    290     ,   290,      290,290,                     246,290      ]
        snums_z1=[290,             290,      290  ,    290     ,   290,      290,290,                     290,290      ]
        snums_z2=[190,             190,      190  ,    190     ,   190,      190,190      ,               204,190]
        snums_z2=[190,             190,      190  ,    190     ,   190,      190,190      ,               190,190]
        snums_z4=[90,             90,      90  ,    90     ,       90,      90,90      ,                   120,90]
        snums_z4=[90,             90,      90  ,    90     ,       90,      90,90      ,                    90,90]
        #snums_z6=[50,             50,      50  ,    50     ,       54,      31,50                               ] #28 m13 problem? (b/f fixed density)
        snums_z6=[50,             50,      50  ,    50     ,       50,      50,50                               ] #28 m13 problem? (b/f fixed density)
        snums_z6=[50,             50,      50  ,    50     ,       50,      50,50            , 50,50                   ] #28 m13 problem? (b/f fixed density)
        #snums_z6=[50,             50,      50  ,    50     ,       54,      50,50                               ] #28 m13 problem? (b/f fixed density)

        if (z_redshift==0.0): snums=snums_z0
        if (z_redshift==0.5): snums=snums_z05
        if (z_redshift==1.0): snums=snums_z1
        if (z_redshift==2.0): snums=snums_z2
        if (z_redshift==4.0): snums=snums_z4
        if (z_redshift==6.0): snums=snums_z6

        if (SHOW_MGAL_SCATTER):
            all_ms=np.zeros((0))
            all_mh=np.zeros((0))

        if (SHOW_MGAL_MHALO_NUMERICS==1):
            if (z_redshift==0.00):
                ## test SPH
                sdirs=['hires_gnew','hires_bh']
                m_primary_list=[12.,    12.,            9.,     9.      ]
                contam_cuts=[0.1,       0.1 ,           1.0,     1.1     ]
                snums=[90,             90,             155,     132     ]
                ## resolution tests
                #sdirs=['m11_gas_new','m11_gas_hr']
                #m_primary_list=[11.,    11.]
                #contam_cuts=[0.5,       0.1 ]
                #snums=[90,             90]
                #snums=[440,            440]
                #sdirs=['m12_gas_lr','m12_mr'] # add new m12_hr as well!
                #m_primary_list=[12.,    12.]
                #contam_cuts=[1.,       0.1 ]
                #snums=[90,             90]
                #sdirs=['m13_gas_lr','m13_gas_mr']
                #m_primary_list=[13.,    13.]
                #contam_cuts=[1.,       0.1 ]
                #snums=[38,             290]
                ## jose's runs, and some of mine, have other numerical variations, add those
                sdirs=['hires_gnew', 'm11_gas_new', 'm12_gas_lr', 'm13_gas_mr', 'hires_bh']
                m_primary_list=[ 12.,   11.,           12.,         13.        ,    12.5] 
                contam_cuts=[ 0.1,     0.5,             1.,         1.         ,     0.01 ]
                snums=[     90,        440,            90,           290       ,    340  ]
                sdirlbl=['m12v (traditional SPH)','m11 (low-res)','m12 (low-res)','m13 (low-res)','m12v (mod. Art. Visc.)']

                sdirs=['m11_gas_new', 'm11_hr_Dec9_2013_','m10_lr_Dec9_2013_','m10_hr_Dec9_2013_']
                m_primary_list=[  12.8,  12.8,                10.,                10.]
                contam_cuts=[ 0.5,      0.1,                0.1,                0.1]
                snums=[    440,         440,                440,                440]
                sdirlbl=[r'm11 (low-res)',r'm11 (standard)',r'm10 (low-res)',r'm10 (standard)']

                sdirs=['m11_gas_new', '','m10_lr_Dec9_2013_','m10_hr_Dec9_2013_']
                sdirs=['m11_hr_Dec9_2013_', 'm11_hhr_Jan9_2013_','m10_lr_Dec9_2013_','m10_hr_Dec9_2013_']
                m_primary_list=[  12.8,  12.8,                10.,                10.]
                m_primary_list=[  13.2,  13.2,                10.,                10.]
                contam_cuts=[ 0.1,      0.1,                0.1,                0.1]
                snums=[    440,         440,                440,                440]
                sdirlbl=[r'm11 (low-res)',r'm11 (standard)',r'm10 (low-res)',r'm10 (standard)']


            if (z_redshift==0.05):
                ## this is the feedback comparison (limited now, add m12/hires_bh with no fb to z=0) ##
                sdirs=['hires_bh_nosn','m10_norad','m11_nofb','z10n192_B1hr']#,'z10n192_B1hr_cwef2en50']
                m_primary_list=[13.,    10.,            11.,       13.     ,            9. ]
                contam_cuts=[0.1,       0.1 ,           1.0,     0.05   ,           0.1     ]
                snums=[124,             90,             155,     441    ,           441     ]
                sdirlbl=['m12v: no SNe','m09: no Rad.','m11: no FB','m12v: no FB']


        f_baryon = 0.162
        marker_syms=['o','s','p','D','*','x','x','D','D']
        colors_vec=['k','b','r','g','c','m','m','y','y']
        marker_syms=['o','s','p','^','*','x','D','+','h','1','v','2','<','3','>','4','H']
        colors_vec=['black','blue','red','green','deepskyblue','darkviolet','orange','magenta','gold','sienna','pink','forestgreen','darkgray']

        if (SHOW_SFR_MAINSEQUENCE):
            i_z_plot=0
            if (z_redshift<=0.5): 
                i_z_plot=1
                zlbl=r'$z<0.5$'
            if (z_redshift==1.0): 
                i_z_plot=2
                zlbl=r'$z\sim 1$'
            if (z_redshift==2.0): 
                i_z_plot=3
                zlbl=r'$z\sim 2$'
            if (z_redshift>=4.0): 
                i_z_plot=4
                zlbl=r'$z\sim 4$'
            plot.subplot(2,4,i_z_plot)
            pylab.axis([7.,12.,1.e-4,1.e3])
            pylab.axis([6.8,11.95,-3.7,2.7])
            pylab.yscale('linear'); pylab.xscale('linear');
            pylab.xlabel(r'$\log(\ M_{\ast}\ /\ M_{\odot}\ )$')
            if (i_z_plot==1):
                pylab.ylabel(r'$\log(\dot{M}_{\ast})\ \  [{\rm M_{\odot}\,yr^{-1}}]$');
            else:
                pylab.ylabel(r' ');
            if (z_redshift != 0) & (z_redshift != 6.0):
                x0_zlvl=10.85
                if(z_redshift==0.5): x0_zlvl=10.5
                pylab.text(x0_zlvl,-3.2,zlbl,fontsize=19,color='k')
            pylab.yticks([2.,1.,0.,-1.,-2.,-3.],['2','1','0','-1','-2','-3'])
            ## z = 0
            if (z_redshift<=0.5): ssfr=np.array([7.e-11,5.e-11,3.e-11]); 
            if (z_redshift==1.0): ssfr=np.array([0.7*9.6e-10,0.9*8.e-10,0.9*4.7e-10]); 
            if (z_redshift==2.0): ssfr=np.array([2.1e-9,1.8e-9,1.4e-9]); 
            if (z_redshift>=4.0): ssfr=np.array([4.5e-9,4.4e-9,4.0e-9]); 
            x=np.arange(8.,10.5,0.1)
            x0=np.array([9.5,10.,10.5])
            y=np.interp(x,x0,np.log10(ssfr)) + x
            pylab.plot(x,y,'k--',linewidth=2.,label=r'Observed')
            if((i_z_plot==1)&(z_redshift==0)):
                pylab.legend(loc='upper left',fontsize=12,frameon=False,labelspacing=0.1)
            ssfr=np.array([7.e-11,5.e-11,3.e-11]); 
            x=np.arange(8.,10.5,0.1)
            x0=np.array([9.5,10.,10.5])
            y=np.interp(x,x0,np.log10(ssfr)) + x
            #pylab.plot(x,y,'k:',linewidth=2.)
            
            if ((z_redshift==2.0) & (2==2)):
                # erb et al.
                x0=[9.57,9.89,10.12,10.31,10.74]
                y0=[1.42,1.60, 1.53, 1.78, 1.96]
                #pylab.plot(x0,y0,marker='o',color='k',linestyle='-')

                # zahid 2012
                x0=[9.6,10.75]; y0=[1.17,2.00];
                x=np.arange(x0[0],x0[1],0.05)
                y=np.interp(x,x0,y0)
                #pylab.plot(x,y,'b--',linewidth=2.)

                # daddi 2007
                x0=[9.26,10.58]; y0=[0.72,1.91];
                x=np.arange(x0[0],x0[1],0.05)
                y=np.interp(x,x0,y0)
                pylab.plot(x,y,'b--',linewidth=2.)

            if ((z_redshift==1.0) & (2==2)):
                # elbaz 2007
                x0=[8.78,10.76]; y0=[-0.34,1.44];
                x=np.arange(x0[0],x0[1],0.05)
                y=np.interp(x,x0,y0)
                pylab.plot(x,y,'b--',linewidth=2.)

            if ((z_redshift==0.0) & (2==2)):
                # sdss
                x0=[9.16,10.43]; y0=[-0.39,0.58];
                x=np.arange(x0[0],x0[1],0.05)
                y=np.interp(x,x0,y0)
                pylab.plot(x,y,'b--',linewidth=2.)

        if (SHOW_MGAL_METALLICITY):
            pylab.axis([3.,12.,-3.5,1.])
            pylab.xscale('linear'); pylab.yscale('linear')

        if (SHOW_KSLAW):
            pylab.axis([-0.7,4.,-5.9,2.])
            pylab.xscale('linear'); pylab.yscale('linear')

        if (SHOW_SSFR):
            msize=2.
            pylab.subplot(2,4,4+2)
            x=np.loadtxt('./ssfr_behroozi_m9.dat')
            #pylab.plot(x[:,0]-1.,x[:,1]+9.,color='k',\
            #    marker='D',markersize=msize,linestyle='',markerfacecolor='none')
            pylab.errorbar(x[:,0]-1.,x[:,1]+9.,yerr=x[:,3],\
                color='k',capsize=0,linestyle='')

            pylab.subplot(2,4,4+3)
            x=np.loadtxt('./ssfr_behroozi_m10.dat')
            #pylab.plot(x[:,0]-1.,x[:,1]+9.,color='k',\
            #    marker='D',markersize=msize,linestyle='',markerfacecolor='none')
            pylab.errorbar(x[:,0]-1.,x[:,1]+9.,yerr=x[:,3],\
                color='k',capsize=0,linestyle='')
            pylab.subplot(2,4,4+3)
            x=np.loadtxt('./ssfr_behroozi_m105.dat')
            #pylab.plot(x[:,0]-1.,x[:,1]+9.,color='k',\
            #    marker='D',markersize=msize,linestyle='',markerfacecolor='none')
            pylab.errorbar(x[:,0]-1.,x[:,1]+9.,yerr=x[:,3],\
                color='k',capsize=0,linestyle='')

            pylab.subplot(2,4,4+4)
            pylab.xlabel(r'${\rm Redshift}\ z$')
            pylab.ylabel(r' ');
            x=np.loadtxt('./ssfr_behroozi_m105.dat')
            #pylab.plot(x[:,0]-1.,x[:,1]+9.,color='k',\
            #    marker='D',markersize=msize,linestyle='',markerfacecolor='none')
            pylab.errorbar(x[:,0]-1.,x[:,1]+9.,yerr=x[:,3],\
                color='k',capsize=0,linestyle='',label=r'Observed')
            if (z_redshift==0.): 
                pylab.legend(loc='lower right',fontsize=12,frameon=False,
                    labelspacing=0.1,numpoints=1)

            pylab.subplot(2,4,4+3)
            pylab.xlabel(r'${\rm Redshift}\ z$')
            pylab.ylabel(r' ');
            x=np.loadtxt('./ssfr_time.dat')
            pylab.scatter((10.**x[:,0])-1.,x[:,1],color='k',\
                marker='',facecolor='none')
            pylab.errorbar((10.**x[:,0])-1.,x[:,1],yerr=0.*x[:,1]+0.35,\
                color='k',capsize=0,linestyle='')

            pylab.subplot(2,4,4+2)
            pylab.xlabel(r'${\rm Redshift}\ z$')
            pylab.ylabel(r' ');
            x=np.loadtxt('./ssfr_time.dat')
            pylab.scatter((10.**x[:,0])-1.,x[:,1],color='k',\
                marker='',facecolor='none')
            pylab.errorbar((10.**x[:,0])-1.,x[:,1],yerr=0.*x[:,1]+0.35,\
                color='k',capsize=0,linestyle='')

            pylab.subplot(2,4,4+1)
            pylab.xlabel(r'${\rm Redshift}\ z$')
            pylab.ylabel(r'$\log( \dot{M}_{\ast}/M_{\ast} )\ \ [{\rm Gyr}^{-1}]$');
            lbl_vec_tmp=[r'$\log(M_{\ast})\sim 8$',r'$\log(M_{\ast})\sim 9$',
                r'$\log(M_{\ast})\sim 10$',r'$\log(M_{\ast})\sim 11$']

            for i_tmp,lbl_tmp in zip([1,2,3,4],lbl_vec_tmp):
                pylab.subplot(2,4,4+i_tmp)
                pylab.yticks([2.,1.,0.,-1.,-2.],['2','1','0','-1','-2'])
                if (z_redshift==0.): 
                    pylab.text(0.,1.4,lbl_tmp,color='k',fontsize=16)


        if (SHOW_MGAL_MHALO_ALLZ):
            #plot.subplot(2,6,i_z+1)
            plot.subplot(2,3,i_z+1)
            all_ms=np.zeros((0))
            all_mh=np.zeros((0))
            #pylab.axis([6.9,14.2,3.,12.5]); pylab.yscale('linear');
            pylab.axis([6.9,13.3,3.3,11.5]); pylab.yscale('linear'); #axis for sim set 1
            pylab.axis([7.9,13.3,3.3,11.8]); pylab.yscale('linear');
            if ((i_z==0)|(i_z==3)):
                pylab.ylabel(r'$\log(\ {\rm M}_{\ast}\ /\ {\rm M}_{\odot}\ )$',fontsize=charsi)
            else:
                pylab.ylabel(' ')
            if ((i_z>=3)):
                pylab.xlabel(r'$\log(\ {\rm M}_{\rm halo}\ /\ {\rm M}_{\odot}\ )$',fontsize=charsi)
            else:
                pylab.xlabel(' ')
            x=np.arange(5.,15.,0.1); y=x+np.log10(0.162) ## universal baryon fraction
            pylab.plot(x,y,'k:')
            ## behroozi relation
            fname='./sm_mh_rel/c_smmr_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'
            v = np.loadtxt(fname)
            x0=np.arange(5.,15.,0.1)
            x=v[:,0]; y=x+v[:,1]; ym=y-v[:,3]; yp=y+v[:,2]
            #pylab.plot(x,y,'m:')
            ## behroozi at the z of interest:
            if (z_redshift==0.0): zkey='0.10'
            if (z_redshift==0.5): zkey='0.10'
            if (z_redshift==1.0): zkey='1.00'
            if (z_redshift==2.0): zkey='2.00'
            if (z_redshift==4.0): zkey='4.00'
            if (z_redshift==6.0): zkey='6.00'
            fname='./sm_mh_rel/c_smmr_z'+zkey+'_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'
            v = np.loadtxt(fname)
            x0=np.arange(5.,15.,0.1)
            x=v[:,0]; y=x+v[:,1]; ym=y-v[:,3]; yp=y+v[:,2]
            #pylab.plot(x,y,'m-')
            pylab.plot(x,yp+0.15,'c')
            pylab.plot(x,ym-0.15,'c')

            ## moster at the z of interest
            zx=z_redshift/(1.+z_redshift)
            m10=11.590; m11=1.195; n10=0.0351; n11=-0.0247; b10=1.376; b11=-0.826; g10=0.608; g11=0.329;
            dm10=0.236; dm11=0.353; dn10=0.0058; dn11=0.0069; db10=0.153; db11=0.225; dg10=0.059; dg11=0.173;
            gz=g10+g11*zx; bz=b10+b11*zx; nz=n10+n11*zx; mz=m10+m11*zx;
            M=np.arange(5.,15.,0.1); x=10.**(M-mz)
            msmin = 0.
            if (z_redshift==0.): msmin=7.5
            if (z_redshift==0.5): msmin=8.2
            if (z_redshift==1.): msmin=8.5
            if (z_redshift==2.): msmin=9.2
            if (z_redshift==4.): msmin=10.
            if (z_redshift==6.): msmin=10.
            ms = M + np.log10(2.*nz) - np.log10(x**gz + x**(-bz))
            #pylab.plot(M,ms,linewidth=2.,color='c',linestyle='-')

            gz=g10+g11*zx; bz=(b10-db10)+(b11-db11)*zx; nz=(n10+dn10)+(n11+dn11)*zx; mz=(m10-dm10)+(m11-dm11)*zx;
            M=np.arange(5.,15.,0.05); x=10.**(M-mz)
            ms = M + np.log10(2.*nz) - np.log10(x**gz + x**(-bz))
            ok=(ms>msmin)
            pylab.plot(M[ok],ms[ok]+0.15,linewidth=1.,color='m',linestyle='-')

            gz=g10+g11*zx; bz=(b10+db10)+(b11+db11)*zx; nz=(n10-dn10)+(n11-dn11)*zx; mz=(m10+dm10)+(m11+dm11)*zx;
            M=np.arange(5.,15.,0.05); x=10.**(M-mz)
            ms = M + np.log10(2.*nz) - np.log10(x**gz + x**(-bz))
            ok=(ms>msmin)
            pylab.plot(M[ok],ms[ok]+0.15,linewidth=1.,color='m',linestyle='-')

            if (z_redshift==0.): zlab=r'$z=0$'
            if (z_redshift==0.5): zlab=r'$z=0.5$'
            if (z_redshift==1.): zlab=r'$z=1$'
            if (z_redshift==2.): zlab=r'$z=2$'
            if (z_redshift==4.): zlab=r'$z=4$'
            if (z_redshift==6.): zlab=r'$z=6$'
            #pylab.text(12.5,3.4,zlab,color='k',fontsize=18)
            #pylab.text(11.65,3.6,zlab,color='k',fontsize=18) #old sims
            pylab.text(11.95,3.6,zlab,color='k',fontsize=18)

            if(2==0):
                plot.subplot(2,6,6+i_z+1)
                pylab.xlabel(r'$\log(\ {\rm M}_{\ast}\ /\ {\rm M}_{\odot}\ )$',fontsize=charsi)
                if (i_z==0):
                    pylab.ylabel(r'$\log(\ {\rm d}N/{\rm d}\log{M_{\ast}}\ \ [{\rm Mpc}^{-3}])$',fontsize=charsi)
                else:
                    pylab.ylabel(' ')
            


        if (SHOW_MGAL_MHALO_Z0):
            plot.figure(1,figsize=(8.,12.))
        
            charsi=20.
            matplotlib.rcParams.update({'font.size':18})

            plot.subplot(211)
            #pylab.axis([6.9,14.2,3.,12.5]); pylab.yscale('linear');
            #pylab.axis([7.7,13.3,3.3,11.5]); pylab.yscale('linear');
            pylab.axis([7.9,13.3,3.1,11.9]); pylab.yscale('linear');

            #pylab.axis([7.9,13.7,3.1,12.4]); pylab.yscale('linear');

            pylab.xlabel(' ')
            pylab.ylabel(r'$\log(\ {\rm M}_{\ast}\ /\ {\rm M}_{\odot}\ )$',fontsize=charsi)
            x=np.arange(5.,15.,0.1); y=x+np.log10(0.162) ## universal baryon fraction
            pylab.plot(x,y,'k:')
            
            #x=np.arange(9.6,11.1,0.05)
            #pylab.plot(x,7.42+(9.28-7.42)/(11.45-9.58)*(x-9.58),'r--')
            #pylab.plot(x,5.85+(9.28-5.65)/(11.45-9.79)*(x-9.79),'r--')

            plot.subplot(212)
            #pylab.axis([6.9,14.2,-4,0.1]); pylab.yscale('linear');
            #pylab.axis([7.7,13.3,-4.0,0.1]); pylab.yscale('linear');
            pylab.axis([7.9,13.3,-4.3,0.1]); pylab.yscale('linear');

            #pylab.axis([7.9,13.7,-4.3,0.1]); pylab.yscale('linear');

            pylab.xlabel(r'$\log(\ {\rm M}_{\rm halo}\ /\ {\rm M}_{\odot}\ )$',fontsize=charsi)
            pylab.ylabel(r'$\log(\ {\rm M}_{\ast}\ /\ f_{b}\,{\rm M}_{\rm halo}\ )$',fontsize=charsi)
            pylab.plot(x,0.*x+0.,'k:',label=r'$M_{\ast}=f_{b}\,{\rm M}_{\rm halo}$')

            ## behroozi relation
            fname='./sm_mh_rel/c_smmr_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'
            v = np.loadtxt(fname)
            x0=np.arange(5.,15.,0.1)
            x=v[:,0]; y=x+v[:,1]; ym=y-v[:,3]; yp=y+v[:,2]
            #pylab.plot(x,y,'m')
            plot.subplot(211)
            pylab.plot(x,yp,'c')
            pylab.plot(x,ym,'c')
            plot.subplot(212)
            pylab.plot(x,(yp-x-np.log10(f_baryon)),'c')
            pylab.plot(x,(ym-x-np.log10(f_baryon)),'c',label=r'Observed (Behroozi)')
            ## interpolate beyond
            xmin = x[0]
            y=yp; x0=np.arange(5.,x[1],0.1); yyp=y[0]+(y[1]-y[0])/(x[1]-x[0])*(x0-x[0])
            y=ym; x0=np.arange(5.,x[1],0.1); yym=y[0]+(y[1]-y[0])/(x[1]-x[0])*(x0-x[0])
            plot.subplot(211)
            pylab.plot(x0,yyp,'c--')
            pylab.plot(x0,yym,'c--')
            plot.subplot(212)
            pylab.plot(x0,(yyp-x0-np.log10(f_baryon)),'c--')
            pylab.plot(x0,(yym-x0-np.log10(f_baryon)),'c--')

            ## moster relation
            x=np.arange(5.,15.,0.1);
            p = [0.02891, 11.866, 1.110, 0.648]
            p = [0.02820, 11.884, 1.057, 0.556]
            p = [0.02820+0.00065, 11.884-0.030, 1.057-0.054, 0.556+0.017]
            m_M = 2.*p[0] / (10.**((x-p[1])*(-p[2])) + 10.**((x-p[1])*(+p[3])))
            y = x + np.log10(m_M)
            ok = (x > xmin)
            plot.subplot(211)
            pylab.plot(x[ok],y[ok]+0.15,'m-')
            pylab.plot(x[ok],y[ok]-0.15,'m-')
            plot.subplot(212)
            pylab.plot(x[ok],(y[ok]+0.15-x[ok]-np.log10(f_baryon)),'m-')
            pylab.plot(x[ok],(y[ok]-0.15-x[ok]-np.log10(f_baryon)),'m-',label=r'Observed (Moster)')
            ok = (x <= xmin)
            plot.subplot(211)
            pylab.plot(x[ok],y[ok]+0.15,'m--')
            pylab.plot(x[ok],y[ok]-0.15,'m--')
            plot.subplot(212)
            pylab.plot(x[ok],(y[ok]+0.15-x[ok]-np.log10(f_baryon)),'m--')
            pylab.plot(x[ok],(y[ok]-0.15-x[ok]-np.log10(f_baryon)),'m--')

            pylab.legend(loc='lower right',fontsize=14,frameon=False,labelspacing=0.2)

        
        if (SHOW_MGAL_MHALO_NUMERICS):
            plot.figure(1,figsize=(8.,12.))
            charsi=20.
            matplotlib.rcParams.update({'font.size':18})
        
            if (z_redshift == 0.00): 
                plot.subplot(2,1,1)
                pylab.xlabel(r' ',fontsize=charsi)
            if (z_redshift == 0.05): 
                plot.subplot(2,1,2)
                pylab.xlabel(r'$\log(\ {\rm M}_{\rm halo}\ /\ {\rm M}_{\odot}\ )$',fontsize=charsi)
            #pylab.axis([7.,14.2,-4.05,0.1])
            #pylab.axis([8.,12.2,-3.2,0.1])
            pylab.axis([9.,12.2,-2.9,0.1])

            f_bar = 0.162
            pylab.yscale('linear');
            pylab.ylabel(r'$\log(\ {\rm M}_{\ast}\ /\ f_{b}\,{\rm M}_{\rm halo}\ )$',fontsize=charsi)

            x=np.arange(5.,15.,0.1); y=x+np.log10(0.162) ## universal baryon fraction
            pylab.plot(x,y*0.,'k:')
            
            ## behroozi relation
            fname='./sm_mh_rel/c_smmr_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'
            v = np.loadtxt(fname)
            x0=np.arange(5.,15.,0.1)
            x=v[:,0]; y=x+v[:,1]; ym=y-v[:,3]; yp=y+v[:,2]
            #pylab.plot(x,y,'m')
            pylab.plot(x,(yp-x-np.log10(f_baryon)),'c')
            pylab.plot(x,(ym-x-np.log10(f_baryon)),'c')
            ## interpolate beyond
            xmin = x[0]
            y=yp; x0=np.arange(5.,x[1],0.1); yyp=y[0]+(y[1]-y[0])/(x[1]-x[0])*(x0-x[0])
            y=ym; x0=np.arange(5.,x[1],0.1); yym=y[0]+(y[1]-y[0])/(x[1]-x[0])*(x0-x[0])
            pylab.plot(x0,(yyp-x0-np.log10(f_baryon)),'c--')
            pylab.plot(x0,(yym-x0-np.log10(f_baryon)),'c--')

            ## moster relation
            x=np.arange(5.,15.,0.1);
            p = [0.02891, 11.866, 1.110, 0.648]
            p = [0.02820, 11.884, 1.057, 0.556]
            p = [0.02820+0.00065, 11.884-0.030, 1.057-0.054, 0.556+0.017]
            m_M = 2.*p[0] / (10.**((x-p[1])*(-p[2])) + 10.**((x-p[1])*(+p[3])))
            y = x + np.log10(m_M)
            ok = (x >= xmin)
            pylab.plot(x[ok],(y[ok]+0.15-x[ok]-np.log10(f_baryon)),'m-')
            pylab.plot(x[ok],(y[ok]-0.15-x[ok]-np.log10(f_baryon)),'m-')
            ok = (x <= xmin+0.1)
            pylab.plot(x[ok],(y[ok]+0.15-x[ok]-np.log10(f_baryon)),'m--')
            pylab.plot(x[ok],(y[ok]-0.15-x[ok]-np.log10(f_baryon)),'m--')
        
            if (z_redshift==0.00):
                ## jose runs, m09
                ms = np.log10(np.array([7.e5, 1.2e6, 1.5e6, 3.4e6, 5.0e6]))
                mh = 0.*ms + np.log10(3.2e9) + 0.05*2.*(np.random.rand(ms.size)-0.5)
                ms_mh = ms-(mh+np.log10(0.162))
                labeler=np.array([r'm09 ($n_{\rm crit}=10$, traditional SPH)',\
                    r'm09 ($n_{\rm crit}=1000$, traditional SPH)',\
                    r'm09 ($n_{\rm crit}=1000$, low-res)',\
                    r'm09 ($n_{\rm crit}=10$, low-res)',\
                    r'm09 (comoving softening, low-res)'])
                for i in range(mh.size):
                    #pylab.plot(mh,ms_mh,marker='D',color='y',linestyle='')
                    marker_sym=marker_syms[i+4]
                    color_sym=colors_vec[4]
                    marksize=10.0
                    pylab.plot([mh[i]],[ms_mh[i]],\
                        marker=marker_sym,markerfacecolor='w',\
                        markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2,\
                        label=labeler[i],linestyle='')

                ms = np.log10(np.array([2.5,1.5,0.8,3.0,1.53,2.6,2.2])*1.e10)
                mh = np.log10(np.array([55.,56.,57.,59.,55.,57.,58.])*0.95*1.e10)
                ms_mh = ms-(mh+np.log10(0.162))
                labeler=np.array([r'm12v (standard)',\
                    r'm12v (traditional SPH)',\
                    r'm12v (mod. RP algorithm)',\
                    r'm12v (mod. RP+SNe algorithm)',\
                    r'm12v (stricter virial SF criterion)',\
                    r'm12v (mod. SNe coupling algorithm)',\
                    r'm12v (mod. Art. Visc.)'])
                for i in range(mh.size):
                    marker_sym=marker_syms[i+0]
                    color_sym=colors_vec[0]
                    marksize=10.0
                    pylab.plot([mh[i]],[ms_mh[i]],\
                        marker=marker_sym,markerfacecolor='w',\
                        markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2,\
                        label=labeler[i],linestyle='')

                ms = np.log10(np.array([8.4,5.4])*1.e10)
                mh = np.log10(np.array([119.,127.])*1.e10)
                ms_mh = ms-(mh+np.log10(0.162))
                labeler=np.array([r'm12i (low-res)',\
                    r'm12i (standard)'])
                for i in range(mh.size):
                    marker_sym=marker_syms[1-i]
                    color_sym=colors_vec[6]
                    marksize=10.0
                    pylab.plot([mh[i]],[ms_mh[i]],\
                        marker=marker_sym,markerfacecolor='w',\
                        markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2,\
                        label=labeler[i],linestyle='')



            if (z_redshift==0.05):
                ## jose runs, m10
                ms = np.log10(np.array([6.e8]))
                mh = 0.*ms + np.log10(np.array(0.99e10)) + 0.05*2.*(np.random.rand(ms.size)-0.5)
                ms_mh = ms-(mh+np.log10(0.162))
                ms_mh = np.array([-2.62,-0.30,-0.43,-2.04,-0.66,-1.93,-1.80])
                mh = np.log10(np.array([0.7e10,0.9e10,0.85e10,0.72e10,0.75e10,0.66e10,0.68e10,0.73e10]))
                labeler=np.array(['All FB on','No FB','No RP+HII','No SW','No SNe','No RP','No HII'])
                labeler=np.array(['All Feedback Active','No Feedback','No Radiation Pressure or HII','No Stellar Wind-Heating','No Supernovae','No Radiation Pressure','No HII (Local Photo-Heating)'])
                for i in range(ms_mh.size):
                    #pylab.plot(mh,ms_mh,marker='D',color='y',linestyle='')
                    marker_sym=marker_syms[i]
                    color_sym=colors_vec[1]
                    marksize=10.0
                    pylab.plot([mh[i]],[ms_mh[i]],\
                        marker=marker_sym,markerfacecolor='w',\
                        markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2,\
                        label=labeler[i],linestyle='')

                ## m12q_noRP (m12q_nov22_crad); 
                ms = np.log10(np.array([9.6e10,2.e9,2.5e10,2.1e10,3.4e10,2.7e10]))
                mh = 0.*ms + np.log10(np.array([0.95e12,1.9e11,5.e11,1.e12,4.e11,4.5e11]))
                ms_mh = ms-(mh+np.log10(0.162))
                colornum = [2, 3, 0, 2 , 0, 0]
                markernum = [5, 0, 0, 0 , 2, 3] 
                for i in range(ms_mh.size):
                    marker_sym=marker_syms[markernum[i]]
                    color_sym=colors_vec[colornum[i]]
                    marksize=10.0
                    pylab.plot([mh[i]],[ms_mh[i]],\
                        marker=marker_sym,markerfacecolor='w',\
                        markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2,\
                        linestyle='')


        if (SHOW_MGAL_MHALO_ALTSIMS):
            SHOW_MGAL_MHALO_ALTSIMS_SHOWOURSIMS = 1
            plot.figure(1,figsize=(8.,12.))
            charsi=20.
            matplotlib.rcParams.update({'font.size':18})

            #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_new', 'm09_dw', 'm13_gas_lr','m13_gas_lr','m14_lr','m14_lr']
            #m_primary_list=[12.,  10.,      12.2,        11.,            9.,     10.,10.,                14.2,    14.9]

            #bad=(sdirs=='hires_bh')
            #m_primary_list[bad]=12.5
            #m_primary_list[bad]=13.

            if (z_redshift<=0.5):
                plot.subplot(2,1,1)
                #pylab.axis([9.,13.,-3.05,0.1])
                if (SHOW_MGAL_MHALO_ALTSIMS_SHOWOURSIMS):
                    pylab.axis([9.,13.,-2.7,0.05])
                    pylab.text(11.05,-2.55,r'$z=0-0.5$',color='k',fontsize=22)
                else:
                    pylab.axis([9.2,13.,-2.3,0.05])
                    pylab.text(11.05,-2.19,r'$z=0-0.5$',color='k',fontsize=22)
                #pylab.text(12.13,-1.83,r'$z=0-0.5$',color='k',fontsize=22)
                #pylab.text(12.13,-2.9,r'$z=0-0.5$',color='k',fontsize=22)
                #pylab.text(10.9,-2.9,r'$z=0-0.5$',color='k',fontsize=22)
                pylab.xlabel(r' ',fontsize=charsi)
            else:
                plot.subplot(2,1,2)
                #pylab.axis([8.3,12.3,-3.05,0.1])
                if (SHOW_MGAL_MHALO_ALTSIMS_SHOWOURSIMS):
                    pylab.axis([8.3,12.15,-2.6,0.05])
                    pylab.text(11.35,-2.5,r'$z=2-3$',color='k',fontsize=22)
                else:
                    pylab.axis([8.3,12.15,-1.9,0.05])
                    pylab.text(8.4,-1.8,r'$z=2-3$',color='k',fontsize=22)
                #m_primary_list[sdirs=='hires_bh']=13.
                #m_primary_list[sdirs=='m11_gas_hr']=15.
                #m_primary_list[sdirs=='m14_lr']=14.2
                #sdirs=['hires_bh',  'zoom_dw', 'm12_mr', 'm11_gas_hr', 'm09_dw', 'm13_gas_lr','m13_gas_mr','m14_lr']#,'m14_lr']
                #m_primary_list=[13.,  10.,      12.2,        12.,            9.,     12.3,12.5,                14.5,    14.9]
                #sdirs=[default_B1,  default_m10, default_m12, default_m11, default_m09, default_m13]#,'m13_gas_lr','m14_lr','m14_lr']
                #m_primary_list=[13.,  10.,      12.2,        12.,            9.,     13.3,12.5,                14.5,    14.9]
                #pylab.text(11.5,-2.9,r'$z=2-3$',color='k',fontsize=22)
                pylab.xlabel(r'$\log(\ {\rm M}_{\rm halo}\ /\ {\rm M}_{\odot}\ )$',fontsize=charsi)
                
            f_bar = 0.162
            pylab.yscale('linear');
            pylab.ylabel(r'$\log(\ {\rm M}_{\ast}\ /\ f_{b}\,{\rm M}_{\rm halo}\ )$',fontsize=charsi)

            x=np.arange(5.,15.,0.1); y=x+np.log10(0.162) ## universal baryon fraction
            pylab.plot(x,y*0.,'k:')
            

            if (z_redshift<=0.5):                
                ## behroozi relation
                fname='./sm_mh_rel/c_smmr_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'
                v = np.loadtxt(fname)
                x0=np.arange(5.,15.,0.1)
                x=v[:,0]; y=x+v[:,1]; ym=y-v[:,3]; yp=y+v[:,2]
                pylab.plot(x,(yp-x-np.log10(f_baryon)),'c')
                pylab.plot(x,(ym-x-np.log10(f_baryon)),'c')
                ## interpolate beyond
                xmin = x[0]
                y=yp; x0=np.arange(5.,x[1],0.1); yyp=y[0]+(y[1]-y[0])/(x[1]-x[0])*(x0-x[0])
                y=ym; x0=np.arange(5.,x[1],0.1); yym=y[0]+(y[1]-y[0])/(x[1]-x[0])*(x0-x[0])
                pylab.plot(x0,(yyp-x0-np.log10(f_baryon)),'c--')
                pylab.plot(x0,(yym-x0-np.log10(f_baryon)),'c--')

                ## moster relation
                x=np.arange(5.,15.,0.1);
                p = [0.02891, 11.866, 1.110, 0.648]
                p = [0.02820, 11.884, 1.057, 0.556]
                p = [0.02820+0.00065, 11.884-0.030 + np.log10(1./0.9), 1.057-0.054, 0.556+0.017]
                m_M = 2.*p[0] / (10.**((x-p[1])*(-p[2])) + 10.**((x-p[1])*(+p[3])))
                y = x + np.log10(m_M)
                ok = (x > xmin)
                pylab.plot(x[ok],(y[ok]+0.125-x[ok]-np.log10(f_baryon)),'m-')
                pylab.plot(x[ok],(y[ok]-0.15-x[ok]-np.log10(f_baryon)),'m-')
                ok = (x <= xmin+0.1)
                pylab.plot(x[ok],(y[ok]+0.125-x[ok]-np.log10(f_baryon)),'m--')
                pylab.plot(x[ok],(y[ok]-0.15-x[ok]-np.log10(f_baryon)),'m--')

                x=[11.90]; y=[10.58] # Guedes 2011
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='orange',linestyle='',label=r'Guedes 2011')
                x=[9.31]; y=[7.24] # Mashchenko 2008
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='m',linestyle='',label=r'Mashchenko 2008')
                x=[9.70,9.94,10.14]; y=[7.90,8.35,8.60] # Stinson 2007
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='c',linestyle='',label=r'Stinson 2007')
                x=[9.61,9.64]; y=[8.69,8.76] # Valcke 2008
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='r',linestyle='',label=r'Valcke 2008')
                x=[10.30,10.54]; y=[8.73,8.71] # Governato 2010
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='burlywood',linestyle='',label=r'Governato 2010')
                x=[9.92,9.94,9.95,9.96,10.01]; y=[7.74,7.80,7.71,7.94,8.02] # Sawala 2011
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='g',linestyle='',label=r'Sawala 2011')
                x=[10.17]; y=[8.29] # Pelupessy 2004
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='k',linestyle='',label=r'Pelupessy 2004')
                f_bar=np.array(0.1622)
                x=[11.13,11.33,11.68,11.71,11.74,11.76]; y=-np.array([0.79,0.74,0.27,0.35,0.37,0.43])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='s',color='m',linestyle='',label=r'Governato 2012') # governato 2012
                x=[11.70,11.76]; y=-np.array([0.14,0.16])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='s',color='c',linestyle='',label=r'Brooks 2011') # brooks 2011
                x=[11.82,11.86,11.88,11.90,12.01,12.12,12.16,12.18]; 
                y=-np.array([0.48,0.60,0.42,0.51,0.36,0.44,0.61,0.46])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='s',color='g',linestyle='',label=r'Scannapieco 2011') # scannapieco 2011
                x=[ 11.864,  11.966,  12.045,  12.058,  12.12 ,  12.162,  12.183, \
                    12.199,  12.265,  12.27 ,  12.344,  12.419,  12.451,  12.503, \
                    12.537,  12.538,  12.632,  12.687,  12.69 ,  12.75 ,  12.809, \
                    12.859,  12.892,  12.953,  12.987,  13.026,  13.026]; 
                y=np.array([-0.3857, -0.4569, -0.5338, -0.3108, -0.4081, -0.7225, -0.5289, \
                    -0.4058, -0.6616, -0.4914, -0.5448, -0.6804, -0.5426, -0.455 , \
                    -0.6343, -0.5961, -0.5909, -0.6471, -0.8614, -0.63  , -0.8036, \
                    -0.8803, -0.6956, -1.016 , -0.7051, -0.8052, -0.6525])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='s',color='burlywood',linestyle='',label=r'Oser 2010') # oser 2010

                x=[ 10.543,  10.65 ,  10.745,  10.864,  10.927,  11.007,  11.101,\
                    11.122,  11.185,  11.276,  11.364,  11.455,  11.574,  11.673,\
                    11.78 ,  11.969,  12.15 ,  12.307,  12.417,  12.488,  12.606,\
                    12.716,  12.822,  12.92 ,  13.082]; 
                y=np.array([-1.787 , -1.724 , -1.636 , -1.564 , -1.476 , -1.412 , -1.384 ,\
                    -1.26  , -1.196 , -1.138 , -1.027 , -0.9169, -0.8268, -0.7541,\
                    -0.6903, -0.7094, -0.7841, -0.8733, -0.9357, -1.015 , -1.086 ,\
                    -1.161 , -1.241 , -1.306 , -1.404 ])+x+np.log10(f_bar)
                #pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='.',color='burlywood',linestyle='-') # guo 2011 (SAM)

                x=np.log10(1./0.7*1.e10*np.array([904.7,657.7,681.3,662.1,640.3,695.4,\
                555.2,546.7,504.4,465.6,462.0,442.4,326.7,381.1,304.6,292.7,261.3,245.8,\
                196.0,189.1,169.4,164.8,161.1,131.4,141.8,141.4,148.3,149.1,124.3,112.0,132.3,93.03,\
                81.68,82.87,56.10,33.22,21.81,40.49,34.79,31.10,28.46,25.74,43.18,22.93,17.51]))
                y=np.log10(1./0.7*1.e10*np.array([36.58,31.92,19.67,26.92,24.56,30.90,\
                17.80,20.65,17.01,20.37,19.98,14.99,15.69,17.46,14.79,10.58,13.37,10.00,\
                10.77,10.93,9.077,8.471,9.379,2.774,6.453,7.375,10.27,8.006,7.461,5.189,7.772,6.250,\
                7.084,2.205,3.710,2.161,0.292,0.732,1.963,1.009,0.579,0.784,0.665,0.373,0.506])) # Hirschmann+2013
                #pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='*',color='blue',markersize=5.,linestyle='',label=r'Hirschmann 2013') # stinson 2013
                x=[12.78,12.22,11.78,11.40,11.62,11.18,10.71]
                y=[0.275,0.401,0.310,0.230,0.272,0.243,0.135]
                pylab.plot(x,np.log10(np.array(y)),marker='*',color='blue',markersize=8.,linestyle='',label=r'Hirschmann 2013') # stinson 2013

                x=np.log10(np.array([39.6,44.5,42.9])*1.e12)
                y=np.log10(np.array([11.0,12.0,15.6])*1.e10)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='*',color='red',markersize=8.,linestyle='',label=r'Feldmann 2010') # stinson 2013

                x=np.array([10.40,10.74,11.09,11.52,11.92,12.30,12.99])
                y=np.array([-2.34,-2.16,-1.99,-1.80,-1.67,-1.62,-2.06])
                pylab.plot(x,np.array(y)-np.log10(f_bar)+0.12,marker='*',color='skyblue',markersize=8.,linestyle='',label=r'De Rossi 2013') # de Rossi 2013

                x=[12.14,12.15]; y=-np.array([0.97,0.80])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='b',linestyle='',markersize=5.,label=r'Okamoto 2012') # okamoto 2012
                x=[11.84]; y=-np.array([0.81])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='r',linestyle='',markersize=5.,label=r'Stinson 2013') # stinson 2013

                x=[11.86,12.02,12.03,12.05,12.20,12.22,12.23]
                y=[10.67,10.63,10.89,10.85,10.81,10.63,11.05] # marinacci+2014
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='yellow',markersize=5.,linestyle='',label=r'Marinacci 2014') # stinson 2013

                x=[11.23,11.48,11.69,11.85,11.92,12.01,12.06,12.11,12.18,12.19,12.21,12.22,12.28,12.33,12.41]
                y=[ 9.53, 9.68,10.00,10.30,10.21,10.66,10.64,10.33,10.73,10.61,11.05,10.56,10.95,10.89,11.29] # aumer+2014
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='lightgreen',markersize=5.,linestyle='',label=r'Aumer 2013') # stinson 2013

                x=np.log10(np.array([3.2e10,2.8e10,2.7e10,2.9e10]))
                y=np.log10(np.array([ 1.3e8, 2.1e7, 7.5e6, 1.6e7])) # trujillo-gomez+2013
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='pink',markersize=5.,linestyle='',label=r'Trujillo-Gomez 2013') # stinson 2013

                x=[11.0]
                y=[18.93] # ceverino+2013
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='magenta',markersize=5.,linestyle='',label=r'Ceverino 2013') # stinson 2013



                #pylab.legend(loc='lower right',fontsize=9,frameon=False,labelspacing=0.1,numpoints=1)
                #pylab.legend(bbox_to_anchor=(1.0,0.515),fontsize=11,frameon=False,labelspacing=0.1,numpoints=1)
                #pylab.legend(bbox_to_anchor=(1.0,0.45),fontsize=11,frameon=False,labelspacing=0.1,numpoints=1)
                pylab.legend(bbox_to_anchor=(1.0,0.62),fontsize=11,frameon=False,labelspacing=0.1,numpoints=1)

            if (z_redshift>=2.0):
                ## behroozi relation
                fname='./sm_mh_rel/c_smmr_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'
                v = np.loadtxt(fname)
                x0=np.arange(5.,15.,0.1)
                x=v[:,0]; y=x+v[:,1]; ym=y-v[:,3]; yp=y+v[:,2]
                pylab.plot(x,y,'m:')
                ## behroozi at the z of interest:
                if (z_redshift==0.0): zkey='0.10'
                if (z_redshift==0.5): zkey='0.10'
                if (z_redshift==1.0): zkey='1.00'
                if (z_redshift==2.0): zkey='2.00'
                if (z_redshift==4.0): zkey='4.00'
                if (z_redshift==6.0): zkey='6.00'
                fname='./sm_mh_rel/c_smmr_z'+zkey+'_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat'
                v = np.loadtxt(fname)
                x0=np.arange(5.,15.,0.1)
                x=v[:,0]; y=x+v[:,1]; ym=y-v[:,3]; yp=y+v[:,2]
                #pylab.plot(x,y-x-np.log10(f_bar),'m-')
                pylab.plot(x-0.1,y+0.15-x-np.log10(f_bar),'c-')
                pylab.plot(x-0.1,y-0.15-x-np.log10(f_bar),'c-')
                xmin = x[0]
                y=yp; x0=np.arange(5.,x[1],0.1); yyp=y[0]+(y[1]-y[0])/(x[1]-x[0])*(x0-x[0])
                y=ym; x0=np.arange(5.,x[1],0.1); yym=y[0]+(y[1]-y[0])/(x[1]-x[0])*(x0-x[0])
                #pylab.plot(x0,(yyp-x0-np.log10(f_bar)),'c--')
                #pylab.plot(x0,(yym-x0-np.log10(f_bar)),'c--')

                ## moster at the z of interest
                zx=z_redshift/(1.+z_redshift)
                m10=11.590; m11=1.195; n10=0.0351; n11=-0.0247; b10=1.376; b11=-0.826; g10=0.608; g11=0.329;
                dm10=0.236; dm11=0.353; dn10=0.0058; dn11=0.0069; db10=0.153; db11=0.225; dg10=0.059; dg11=0.173;
                gz=g10+g11*zx; bz=b10+b11*zx; nz=n10+n11*zx; mz=m10+m11*zx;
                M=np.arange(5.,15.,0.1); x=10.**(M-mz)
                msmin = 0.
                if (z_redshift==0.): msmin=7.5
                if (z_redshift==0.5): msmin=8.2
                if (z_redshift==1.): msmin=8.5
                if (z_redshift==2.): msmin=9.2
                if (z_redshift==2.): msmin=8.8
                if (z_redshift==4.): msmin=10.
                if (z_redshift==6.): msmin=10.
                ms = M + np.log10(2.*nz) - np.log10(x**gz + x**(-bz))
                ok=(ms>msmin)
                #pylab.plot(M[ok],ms[ok]-M[ok]-np.log10(f_bar),linewidth=1.,color='m',linestyle='-')
                pylab.plot(M[ok]-0.2,ms[ok]+0.15-M[ok]-np.log10(f_bar),linewidth=1.,color='m',linestyle='-')
                pylab.plot(M[ok]-0.2,ms[ok]-0.15-M[ok]-np.log10(f_bar),linewidth=1.,color='m',linestyle='-')
                ok = (x <= xmin)
                #pylab.plot(M[ok],ms[ok]+0.15-M[ok]-np.log10(f_bar),linewidth=1.,color='m',linestyle='--')
                #pylab.plot(M[ok],ms[ok]-0.15-M[ok]-np.log10(f_bar),linewidth=1.,color='m',linestyle='--')

                gz=g10+g11*zx; bz=(b10-db10)+(b11-db11)*zx; nz=(n10+dn10)+(n11+dn11)*zx; mz=(m10-dm10)+(m11-dm11)*zx;
                M=np.arange(5.,15.,0.05); x=10.**(M-mz)
                ms = M + np.log10(2.*nz) - np.log10(x**gz + x**(-bz))
                ok=(ms>msmin)
                #pylab.plot(M[ok],ms[ok]+0.15-M[ok]-np.log10(f_bar),linewidth=2.,color='c',linestyle='--')
                gz=g10+g11*zx; bz=(b10+db10)+(b11+db11)*zx; nz=(n10-dn10)+(n11-dn11)*zx; mz=(m10+dm10)+(m11+dm11)*zx;
                M=np.arange(5.,15.,0.05); x=10.**(M-mz)
                ms = M + np.log10(2.*nz) - np.log10(x**gz + x**(-bz))
                ok=(ms>msmin)
                #pylab.plot(M[ok],ms[ok]+0.15-M[ok]-np.log10(f_bar),linewidth=2.,color='c',linestyle='--')

                x=[11.20,10.85]; y=[10.27,9.98] # Guedes 2011
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='orange',linestyle='')
                #x=[9.31]; y=[7.24] # Mashchenko 2008 (don't quote high-z numbers)
                #pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='m',linestyle='')
                #x=[9.70,9.94,10.14]; y=[7.90,8.35,8.60] # Stinson 2007 (don't quote high-z numbers)
                #pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='c',linestyle='')
                x=[9.61,9.64]; y=[8.69,8.76] # Valcke 2008 (assume non-evolving parameters)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='r',linestyle='')
                x=[10.30,10.54]; y=[8.73,8.71] # Governato 2010 (no high-z properties quoted)
                #pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='burlywood',linestyle='')
                x=[8.57,8.48,8.59,8.47]; y=[6.62,6.62,6.83,6.77] # Sawala 2011
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='g',linestyle='')
                #x=[10.17]; y=[8.29] # Pelupessy 2004 (idealized sims, no high-z version)
                #pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='o',color='k',linestyle='')
                f_bar=np.array(0.1622)
                x=[10.57,10.89,10.90,11.21,11.40,11.46,10.26,10.57,10.60,10.70,11.11]; 
                y=-np.array([0.90,0.62,0.99,0.87,1.05,0.57,0.49,0.55,0.65,0.66,0.62])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='s',color='m',linestyle='') # governato 2012
                x=[11.34,11.49,10.64,11.15]; y=-np.array([0.43,0.17,0.19,0.36])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='s',color='c',linestyle='') # brooks 2011
                x=[ 11.497,  11.505,  11.563,  11.651,  11.85 ,  11.917,  11.965,\
                    12.028,  11.019,  11.362,  11.374,  11.468,  11.661,  11.811]
                y=np.array([-0.3258, -0.2818, -0.8017, -0.2739, -0.6302, -0.5484, -0.3696,\
                   -0.37  , -0.909 , -0.3959, -0.7241, -0.2268, -0.7174, -0.3172])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='s',color='g',linestyle='') # scannapieco 2011
                x=[ 11.36 ,  11.454,  11.47 ,  11.567,  11.591,  11.592,  11.635,\
                    11.768,  11.812,  11.867,  11.882,  11.957,  12.023,  12.023,\
                    12.055,  12.074,  12.114,  12.138,  12.161,  12.31 ,  12.33 ,\
                    12.385,  12.455,  12.483,  12.494,  12.558,  12.675,\
                    10.759,  11.098,  11.248,  11.259,  11.362,  11.374,  11.46 ,\
                    11.484,  11.547,  11.57 ,  11.602,  11.629,  11.685,  11.704,\
                    11.704,  11.775,  11.846,  11.96 ,  12.086,  12.192,  12.291]        
                y=np.array([-0.3309, -0.261 , -0.4313, -0.6843, -0.4467, -0.3175, -0.3853,\
                   -0.4419, -0.3452, -0.3808, -0.5188, -0.44  , -0.6312, -0.549 ,\
                   -0.5228, -0.6873, -0.5965, -0.444 , -0.544 , -0.5419, -0.6741,\
                   -0.6363, -0.6543, -0.493 , -0.7749, -0.605 , -0.6556,\
                   -0.245 , -0.5551, -0.4214, -0.6559, -0.601 , -0.4634, -0.2765,\
                   -0.3529, -0.5291, -0.4736, -0.3244, -0.2221, -0.8905, -0.4014,\
                   -0.8056, -0.5806, -0.417 , -0.4003, -0.5096, -0.7184, -0.7806])\
                    +x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='s',color='burlywood',linestyle='') # oser 2010
                x=[ 10.462,  10.477,  10.512,  10.527,  10.56 ,  10.58 ,  10.583,\
                    10.596,  10.638,  10.744,  10.744,  10.804,  10.85 ,  10.937,\
                    11.126,  11.201,\
                    10.011,  10.019,  10.02 ,  10.054,  10.058,  10.074,  10.108,\
                    10.117,  10.117,  10.167,  10.185,  10.193,  10.235,  10.259,\
                    10.306,  10.31 ,  10.319,  10.349,  10.382,  10.394,  10.4  ]
                y=np.array([-1.194 , -1.379 , -1.514 , -1.629 , -1.338 , -1.13  , -1.321 ,\
                    -0.9575, -1.574 , -1.34  , -1.295 , -1.043 , -1.59  , -1.25  ,\
                   -1.321 , -1.16,\
                   -1.23232478, -1.30856485, -1.51870073, -1.74714697, -1.51984927,\
                    -1.34998405, -1.23388472, -1.44514757, -1.61457249, -1.23485321,\
                     -1.59192971, -1.48386142, -1.6064248 , -1.13948232, -1.30768196,\
                     -1.23173198, -1.13596382, -1.58586264, -1.35654732, -1.13691517,\
                     -1.57560845])+x+np.log10(f_bar)
                #pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='.',color='chartreuse',linestyle='') # genel 2012

                x=[11.88,11.97,11.69,11.78]; y=-np.array([0.90,0.64,0.89,0.58])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='b',linestyle='',markersize=5.) # okamoto 2012
                x=[11.29,11.55]; y=-np.array([1.39,1.27])+x+np.log10(f_bar)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='r',linestyle='',markersize=5.) # stinson 2013

                x=[11.60,11.66,11.67,11.75,11.88,11.90,11.96]
                y=[10.46,10.30,10.04, 9.99,10.32,10.24,10.39] # marinacci+2014
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='yellow',linestyle='',markersize=5.)

                x=[10.99,11.04,11.10,11.14,11.25,11.27,11.29,11.31,11.49,11.50,11.64,11.66,11.71,11.80,11.85,11.86,11.92,11.93]
                y=[ 9.05, 9.09, 9.25, 9.03, 9.43, 9.07, 9.27, 9.45, 9.20, 9.70, 9.46,10.08,10.09,10.18,10.03,10.15, 9.93,10.25] # aumer+2014
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='lightgreen',linestyle='',markersize=5.)

                x=[11.0]
                y=[8.93] # ceverino+2013
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='magenta',markersize=5.,linestyle='') # stinson 2013
                x=np.log10(np.array([8.3e10, 0.99e10,1.e10,1.01e10,0.98e10]))
                y=np.log10(np.array([ 1.4e9, 8.e6, 3.4e6, 2.5e6, 1.e6])) # trujillo-gomez+2013
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='D',color='pink',markersize=5.,linestyle='') # stinson 2013
                x=[12.18,11.76,11.25,10.72,11.22,10.70,10.26]
                y=[0.140,0.155,0.061,0.034,0.095,0.040,0.022]
                pylab.plot(x,np.log10(np.array(y)),marker='*',color='blue',markersize=5.,linestyle='') # stinson 2013

                x=np.log10(np.array([1.04,1.82,0.61])*1.e12)
                y=np.log10(np.array([17.8,28.2,16.9])*1.e10)
                pylab.plot(x,np.array(y)-np.array(x)-np.log10(f_bar),marker='*',color='red',markersize=5.,linestyle='') # feldmann 2010

                x=np.array([10.39,10.42,10.70,10.73,11.09,11.11,11.50,11.54,11.79,11.80])
                y=np.array([-2.18,-2.27,-2.08,-1.88,-1.81,-1.62,-1.62,-1.51,-1.56,-1.60])
                pylab.plot(x,np.array(y)-np.log10(f_bar)+0.12,marker='*',color='skyblue',markersize=8.,linestyle='',label=r'De Rossi 2013') # de Rossi 2013

                
        for mh_primary,color_sym,contam_cut,marker_sym,snum,sdir,sdir_lbl in \
                zip(m_primary_list,colors_vec,contam_cuts,marker_syms,snums,sdirs,sdirlbl):
            infiname = local_get_halo_id_filename(snapdir_m+sdir,snum,full_halo_results=1,four_char=0)
            print('loading ',infiname)
            infi=h5py.File(infiname,'r')
            print(list(infi)) ## see available arrays
    
            ## load a brick of quantities
            mh=np.log10(1.e-8 + np.array(infi['halo_m'])) + 10.
            ms=np.log10(1.e-8 + np.array(infi['star_m'])) + 10.
            mg=np.log10(1.e-8 + np.array(infi['gas_m'])) + 10.
            contam = np.array(infi['halo_contamination'])
            ns = np.array(infi['star_n'])
            ng = np.array(infi['gas_n'])

            mh_min = mh_primary - 3.0 
            ## scatter in Mgal-Mhalo seems to blow up between 0.001-0.01*M_primary: real?
            ok_s = (ns > 1) & (contam < contam_cut) & (mh > mh_min)
            ok_g = (ng > 1) & (contam < contam_cut) & (mh > mh_min)

            mh += 0.*0.15*(np.random.rand(mh.size)-0.5)
        
            
            if (SHOW_SSFR):
                marksize=10.
                if (marker_sym=='*'): marksize=20.
                sfr = np.array(infi['gas_sfr'])
                x = ms[ok_s]
                y = np.log10(sfr[ok_s])
                y = np.log10(sfr[ok_s]/(10.**ms[ok_s] / 1.e9))
                
                ms_vec = np.array([8., 9., 10., 11.])
                d_ms = 1.0
                for m_s,i_ms in zip(ms_vec,range(ms_vec.size)):
                    ok = (x >= m_s-d_ms) & (x <= m_s+d_ms)
                    if (x[ok].size > 0):
                        plot.subplot(2,4,4+i_ms+1)
                        pylab.axis([-0.2,6.2,-2.,2.])
                        xx = z_redshift + 0.1*2.*(np.random.rand(x[ok].size)-0.5)
                        pylab.scatter(xx,y[ok],\
                          marker=marker_sym,color=color_sym)
                        
                        
            if (SHOW_SFR_MAINSEQUENCE):
                plot.subplot(2,4,i_z_plot)
                sfr = np.array(infi['gas_sfr'])
                marksize=10.
                if (marker_sym=='*'): marksize=20.
                #pylab.plot(ms,sfr,\
                #    marker=marker_sym,markerfacecolor='w',\
                #    markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2)
                pylab.scatter(ms[ok_s],np.log10(sfr[ok_s]),marker=marker_sym,color=color_sym)

            if (SHOW_KSLAW):
                sfr = np.array(infi['gas_sfr'])
                #mg = np.array(infi['
                marksize=10.
                if (marker_sym=='*'): marksize=20.
                pylab.scatter(ms[ok_s],np.log10(sfr[ok_s]),marker=marker_sym,color=color_sym)

        
            if (SHOW_MGAL_MHALO_Z0):
                plot.subplot(211)
                marksize=10.
                if(sdir_lbl=='m11')&(z_redshift<=0.05):
                    mh[ok_s]+=0.3; # halo finder systematically under-estimated mass
                if(sdir_lbl=='m12i')&(z_redshift<=0.05):
                    mh[ok_s]+=0.3; # halo finder systematically under-estimated mass
                if (marker_sym=='*'): marksize=20.
                main = (ms[ok_s] == np.max(ms[ok_s]))
                if (sdir_lbl=='m14') & (snum==275): main=[]
                q, = pylab.plot(mh[ok_s][main],ms[ok_s][main],linestyle='',\
                        marker=marker_sym,markerfacecolor='w',\
                        markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2,\
                        label=sdir_lbl)
                if (sdir_lbl=='m09'): leg_09=q
                if (sdir_lbl=='m10'): leg_10=q
                if (sdir_lbl=='m11'): leg_11=q
                if (sdir_lbl=='m12q'): leg_12q=q
                if (sdir_lbl=='m12v'): leg_12v=q
                if (sdir_lbl=='m12i'): leg_12i=q
                if (sdir_lbl=='m13'): leg_13=q
                if (sdir_lbl=='m14'): leg_14=q
                pylab.scatter(mh[ok_s],ms[ok_s],marker=marker_sym,color=color_sym)

                plot.subplot(212)
                pylab.scatter(mh[ok_s],(ms[ok_s]-(mh[ok_s]+np.log10(f_baryon))),marker=marker_sym,color=color_sym)
                pylab.plot(mh[ok_s][main],ms[ok_s][main]-(mh[ok_s][main]+np.log10(f_baryon)),\
                    marker=marker_sym,markerfacecolor='w',\
                    markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2)


            if (SHOW_MGAL_MHALO_ALTSIMS):
                if (sdir=='m11_hhr_Jan9_2013_')|(sdir=='m11_gas_new')|(sdir=='m11_nofb')|(sdir=='m11_hr_Dec9_2013_'): 
                        mh+=0.25;
                if (z_redshift<=0.5):
                    plot.subplot(2,1,1)
                else:
                    plot.subplot(2,1,2)
                marksize=5.*2.
                if (marker_sym=='*'): marksize=10.*2.
                main = (ms[ok_s] == np.max(ms[ok_s]))
                if (SHOW_MGAL_MHALO_ALTSIMS_SHOWOURSIMS):
                    if(sdir_lbl=='m12i')&(z_redshift<=0.05):
                        mh[ok_s]+=0.3;
                    if(sdir_lbl=='m11')&(z_redshift<=0.05):
                        mh[ok_s]+=0.1;
                    fbary=f_baryon;
                    xcorr=0.;
                    #if(z_redshift>=2.00): fbary *= 2.;
                    if(z_redshift>=2.00): xcorr = 0.2;
                    pylab.scatter(mh[ok_s]+xcorr,(ms[ok_s]-(mh[ok_s]+xcorr+np.log10(fbary))),marker=marker_sym,color='0.5')
                    pylab.plot(mh[ok_s][main]+xcorr,ms[ok_s][main]-(mh[ok_s][main]+xcorr+np.log10(fbary)),\
                        marker=marker_sym,markerfacecolor='w',\
                        markeredgecolor='0.5',markersize=marksize,markeredgewidth=2)


            if (SHOW_MGAL_MHALO_NUMERICS):
                marksize=10.
                tnum=0; mnum=i;
                if (sdir=='hires_gnew'): mh+=0.5;
                if (sdir=='m12_gas_lr'): mh+=0.5;

                if (sdir=='z10n192_B1hr')|(sdir=='hires_gnew')|(sdir=='hires_bh')|(sdir=='hires_bh_nosn'): 
                    tnum=0;
                    if(sdir=='z10n192_B1hr'): mnum=1;
                    if(sdir=='hires_bh_nosn'): mnum=4;
                if (sdir=='m10_norad'):
                    tnum=4; mnum=2;
                if(sdir=='m10_lr_Dec9_2013_')|(sdir=='m10_hr_Dec9_2013_'): 
                    tnum=1; mnum=2;
                    if(sdir=='m10_lr_Dec9_2013_'): mnum=1;
                    if(sdir=='m10_hr_Dec9_2013_'): mnum=0;
                if (sdir=='m12_gas_lr'): 
                    tnum=2; 
                if (sdir=='m11_hhr_Jan9_2013_')|(sdir=='m11_gas_new')|(sdir=='m11_nofb')|(sdir=='m11_hr_Dec9_2013_'): 
                    tnum=3; mnum=1;
                    if (sdir=='m11_hhr_Jan9_2013_'): 
                        mh+=0.3;
                        mnum=0;
                    if (sdir=='m11_hr_Dec9_2013_'): 
                        mh+=0.25;
                if (sdir=='m13_gas_mr'): 
                    tnum=5;


                if (marker_sym=='*'): marksize=20.
                main = (ms[ok_s] == np.max(ms[ok_s]))
                pylab.scatter(mh[ok_s],(ms[ok_s]-(mh[ok_s]+np.log10(f_baryon))),marker=marker_syms[mnum],color=colors_vec[tnum])
                if(z_redshift<=0.0):
                    pylab.plot(mh[ok_s][main],ms[ok_s][main]-(mh[ok_s][main]+np.log10(f_baryon)),\
                        marker=marker_syms[mnum],markerfacecolor='w',\
                        markeredgecolor=colors_vec[tnum],markersize=marksize,markeredgewidth=2,\
                        label=sdir_lbl,linestyle='')
                else:
                    pylab.plot(mh[ok_s][main],ms[ok_s][main]-(mh[ok_s][main]+np.log10(f_baryon)),\
                        marker=marker_syms[mnum],markerfacecolor='w',\
                        markeredgecolor=colors_vec[tnum],markersize=marksize,markeredgewidth=2,\
                        linestyle='')

            if (SHOW_MGAL_MHALO_ALLZ):
                if(sdir_lbl=='m11')&(z_redshift<=10.):
                    mh[ok_s]+=0.3/(1.+z_redshift); # halo finder systematically under-estimated mass
                if(sdir_lbl=='m12i')&(z_redshift<=0.5):
                    mh[ok_s]+=0.3;
                plot.subplot(2,3,i_z+1)
                marksize=10.
                if (marker_sym=='*'): marksize=20.
                main = (ms[ok_s] == np.max(ms[ok_s]))
                if (sdir_lbl=='m14') & (mh_primary==14.2): main=[]
                pylab.plot(mh[ok_s][main],ms[ok_s][main],\
                    marker=marker_sym,markerfacecolor='w',\
                    markeredgecolor=color_sym,markersize=marksize,markeredgewidth=2)
                pylab.scatter(mh[ok_s],ms[ok_s],marker=marker_sym,color=color_sym)
                np.savetxt('z_'+"%.2f" % z_redshift+'_'+sdir_lbl,[mh[ok_s],ms[ok_s]])
                
                all_mh=np.concatenate((all_mh,np.array(mh[ok_s])))
                all_ms=np.concatenate((all_ms,np.array(ms[ok_s])))
            
            if (SHOW_MGAL_SCATTER):
                all_mh=np.concatenate((all_mh,np.array(mh[ok_s])))
                all_ms=np.concatenate((all_ms,np.array(ms[ok_s])))

            """
            pylab.axis([7.,11.3,1.e-3,100.])
            pylab.yscale('log')
            sfr=np.array(infi['gas_sfr'])
            pylab.scatter(ms[ok_s],sfr[ok_s],marker=marker_sym,color=color_sym)
        
            x=np.arange(3.,12.,0.1)
            y=(10.**0.2) * (10.**(x-10.)) # z=0.02-0.2 (Wuyts et al)
            pylab.plot(x,y,'k-')
            y=(10.**1.2) * (10.**(x-10.)) # z=0.5-1.5 (Wuyts et al)
            pylab.plot(x,y,'b--')
            y=(10.**1.52) * (10.**(x-10.)) # z=1.5-2.5 (Wuyts et al)
            pylab.plot(x,y,'r:')
        
            """
        
            """
            pylab.axis([3.,11.,3.,12.])
            ## this is total gas mass, not gas mass *within* the disk radius, desired here
            ##  to compare with some of the observations
            pylab.scatter(ms[ok_s],mg[ok_s],marker=marker_sym,color=color_sym)
            obs=np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/woo08_mstar_mgas.dat')
            pylab.plot(obs[:,0],obs[:,1],'mo')
            obs=np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/geha06_mstar_mgas.dat')
            pylab.plot(obs[:,0],obs[:,1],'go')
            obs=np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/lee06_dwarfs.dat')
            pylab.plot(obs[:,6],obs[:,1]+obs[:,6],'ro')
            """

            """
            rs = np.array(infi['star_r_half'])
            pylab.scatter(ms[ok_s],rs[ok_s],marker=marker_sym,color=color_sym)
        
            x=np.arange(3.,12.,0.1)
            alpha,beta,gamma,x0=0.14,0.39,0.10,np.log10(3.98e10)
            y=gamma * 10.**(alpha*x) * (1.+10.**(x-x0))**(beta-alpha)
            pylab.plot(x,y)
            pylab.plot(x,y*np.exp(0.47),'b--')
            pylab.plot(x,y*np.exp(-0.47),'b--')
            a,b=0.56,3.47e-6
            y=b*10.**(a*x)
            pylab.plot(x,y,'r')
            pylab.plot(x,y*np.exp(0.47),'r--')
            pylab.plot(x,y*np.exp(-0.47),'r--')
            """
        
            ## these data points are really about the enclosed mass 
            ##  inside the effective radius (from the kinematic models) -- 
            ##   so they *include* dark matter
            #obs=np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/wolf10_size_mass.dat')
            #pylab.plot(obs[:,1],obs[:,0],'ro')
            ## (trickier, need to estimate m_dm inside stellar Re, say; can do with gas 
            ##   as well; but looks like agrees well when done)

            """
            rs = np.array(infi['star_r_half'])
            rh = mycosmo.rvir_of_mhalo(mh[ok_s]+np.log10(0.7),1.)
            #pylab.scatter(mh[ok_s],rs[ok_s],marker=marker_sym,color=color_sym)
            #pylab.scatter(0.015*rh,rs[ok_s],marker=marker_sym,color=color_sym)
            #pylab.scatter(0.035*rh,rs[ok_s],marker=marker_sym,color=color_sym) #z~1-2
            pylab.scatter(0.025*rh,rs[ok_s],marker=marker_sym,color=color_sym) #z~4-6
            x=10.**np.arange(-3.,3.,0.1)
            pylab.plot(x,x,'k--')
            pylab.plot(x,x*10.**0.5,'k:')
            pylab.plot(x,x*10.**(-0.5),'k:')
            """

            #"""
        
            z_solar_all=np.array([ 0.02, 0.28, 3.26e-3, 1.32e-3, 8.65e-3, 2.22e-3, 9.31e-4, \
                1.08e-3, 6.44e-4, 1.01e-4, 1.73e-3 ])
            #"""
            if (SHOW_MGAL_METALLICITY): 
            
                jm = 10 # pick which species
                # C (jm=2): strong increase with mass, halos & gals similar
                #  - broadly similar with N
                # O,Ne,Mg similar (same trend, little more offset: smallest halos slightly O,Ne,Mg-depleted)
                # Si,S,Ca similar (little closer between halo & gal)
                # Fe looks similar in gals, different in halos -- substantially higher Fe abundances?
                z_sun = z_solar_all[jm]
        
                SHOW_STELLAR_METALS_HERE = 1

                main = (ms[ok_s] == np.max(ms[ok_s]))
                if (SHOW_STELLAR_METALS_HERE):
                
                    zs = np.log10(np.array(infi['star_z'][:,jm]) / z_sun) ## 'total' metallicity
                    pylab.scatter(ms[ok_s][main],zs[ok_s][main],marker=marker_sym,color='k')
                    zs = np.log10(np.array(infi['star_z_half'][:,jm]) / z_sun) ## 'total' metallicity (inside R_e)
                    pylab.scatter(ms[ok_s][main],zs[ok_s][main],marker=marker_sym,color='r')

                    ## stellar phase -- compare directly with Fe
                    FeH_obs = np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/woo08_FeH.dat')
                    pylab.plot(FeH_obs[:,0],FeH_obs[:,1],'go')

                    ## this is gas-phase, should really be compared separately 
                    ZH_obs = np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/tremonti04_ZMass.dat')
                    normer = -8.69 # (from Tremonti, this is their unit of solar metallicity)
                    normer = -12.-np.log10(z_solar_all[2])
                    normer = -12.-np.log10(0.001574)
                    pylab.plot(ZH_obs[:,0],ZH_obs[:,3]+normer,'g-')
                    pylab.plot(ZH_obs[:,0],ZH_obs[:,5]+normer,'g--')
                    pylab.plot(ZH_obs[:,0],ZH_obs[:,1]+normer,'g--')

                    ZH_obs = np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/lee06_dwarfs.dat')
                    pylab.plot(ZH_obs[:,6],ZH_obs[:,4]+normer,'mo')

                else:

                    zg = np.log10(np.array(infi['gas_z'][:,jm]) / z_sun) ## 'total' metallicity
                    pylab.scatter(ms[ok_s][main],zg[ok_s][main],marker=marker_sym,color='b')
                    zg = np.log10(np.array(infi['gas_z_half'][:,jm]) / z_sun) ## 'total' metallicity (inside R_e)
                    pylab.scatter(ms[ok_s][main],zg[ok_s][main],marker=marker_sym,color='c')
        
                    ## this is gas-phase, should really be compared separately 
                    ZH_obs = np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/tremonti04_ZMass.dat')
                    normer = -8.69 # (from Tremonti, this is their unit of solar metallicity)
                    normer = -12.-np.log10(z_solar_all[2])
                    normer = -12.-np.log10(0.00574)
                    pylab.plot(ZH_obs[:,0],ZH_obs[:,3]+normer,'g-')
                    pylab.plot(ZH_obs[:,0],ZH_obs[:,5]+normer,'g--')
                    pylab.plot(ZH_obs[:,0],ZH_obs[:,1]+normer,'g--')

                    ZH_obs = np.loadtxt('/Users/phopkins/Documents/work/plots/zooms/lee06_dwarfs.dat')
                    pylab.plot(ZH_obs[:,6],ZH_obs[:,4]+normer,'mo')
        
        
            #"""
            """
            ## good alpha-to-iron: Mg/Fe
            j1 = 6 # Mg
            j2 = 10 # Fe
            # quite flat, though drops a solid amount for most massive system
            #   more interesting is that halos have systematically lower Mg/Fe -- 
            #   alpha elements escaped while 1a's come in later? 

            j1 = 2 # C
            j2 = 10 # Fe 
            # C/Fe systematically rises to smaller systems: trapping more Fe? 
            ## halos also lower C/Fe, but offset switches from highest to lowest-mass, 
            ##  (lowest mass have high C/Fe; all the Fe escaping?)
        
            #j1 = 2 # C
            #j2 = 6 # Mg
            # C/Mg pretty flat, halos much higher in C/Mg: so depleted in Mg?

            z_sun = z_solar_all[j1]/z_solar_all[j2]
            zs = np.log10(np.array(infi['star_z'][:,j1])/np.array(infi['star_z'][:,j2]) / z_sun) ## 'total' metallicity
            pylab.scatter(ms[ok_s],zs[ok_s],marker=marker_sym,color='k')
            zs = np.log10(np.array(infi['star_z_half'][:,j1])/np.array(infi['star_z_half'][:,j2]) / z_sun) ## 'total' metallicity
            pylab.scatter(ms[ok_s],zs[ok_s],marker=marker_sym,color='r')
            zs = np.log10(np.array(infi['gas_z'][:,j1])/np.array(infi['gas_z'][:,j2]) / z_sun) ## 'total' metallicity
            pylab.scatter(ms[ok_s],zs[ok_s],marker=marker_sym,color='b')
            zs = np.log10(np.array(infi['gas_z_half'][:,j1])/np.array(infi['gas_z_half'][:,j2]) / z_sun) ## 'total' metallicity
            pylab.scatter(ms[ok_s],zs[ok_s],marker=marker_sym,color='c')
            """
            """
                recall, element 0 is 'total' metal mass, solar = 0.02: 
                the following 10 are in the following solar abundances & order:
                All.SolarAbundances[0]=0.02;    // total metallicity
                All.SolarAbundances[1]=0.28;    // He
                All.SolarAbundances[2]=3.26e-3; // C
                All.SolarAbundances[3]=1.32e-3; // N
                All.SolarAbundances[4]=8.65e-3; // O
                All.SolarAbundances[5]=2.22e-3; // Ne
                All.SolarAbundances[6]=9.31e-4; // Mg
                All.SolarAbundances[7]=1.08e-3; // Si
                All.SolarAbundances[8]=6.44e-4; // S
                All.SolarAbundances[9]=1.01e-4; // Ca
                All.SolarAbundances[10]=1.73e-3; // Fe
            """        
        
    
            """
            ok=(contam<0.03)
            pylab.scatter(mh[ok],ms[ok],marker=marker_sym,color='k')
            ok=(contam>=0.03) & (contam<0.1)
            pylab.scatter(mh[ok],ms[ok],marker=marker_sym,color='b')
            ok=(contam>=0.1) & (contam<=0.25)
            pylab.scatter(mh[ok],ms[ok],marker=marker_sym,color='g')
            ok=(contam>=0.25) & (contam<=0.5)
            pylab.scatter(mh[ok],ms[ok],marker=marker_sym,color='m')
            ok=(contam>0.5)
            pylab.scatter(mh[ok],ms[ok],marker=marker_sym,color='r')
            """
        
            infi.close()

        if ((SHOW_MGAL_MHALO_ALLZ)&(2==0)):
            ## have all_ms and all_mh
            ok=(all_ms > 2.); ms=all_ms[ok]; mh=all_mh[ok]; 
            s=np.argsort(mh); mh=mh[s]; ms=ms[s];
            dmh=0.05; dms=1.*dmh; 
            mh_grid = np.arange(mh[0],mh[-1],dmh); 
            
            popt,pcov = optimize.curve_fit(fitfun_poly5,mh,ms,p0=[1.,1.,1.,1.,1.,1.]);
            #pylab.plot(mh_grid,fitfun_poly5(mh_grid,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),'r-')
            ms_expected = fitfun_poly5(mh,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])       
            
            ms_std = 0.*mh;
            d_i = 5
            d_i = 2
            for i in range(mh.size):
                imin=i-d_i; 
                if (imin<=0): imin=0
                imax=i+d_i; 
                if (imax>=mh.size): imax=mh.size
                d_ms_in=ms[imin:imax] - ms_expected[imin:imax]
                ms_std[i] = np.std(d_ms_in)/np.sqrt(2.)

            plot.subplot(2,6,i_z+1+6)
            ms_grid = np.arange(3.,13.,dms);
            zz = z_redshift
            if (zz==6.0): zz+=2.0
            num_mhalo = mycosmo.mass_function( mh_grid+np.log10(0.7) , zz )
            num_ms = 0.*ms_grid
            ms_std_grid = np.interp(mh_grid,mh,ms_std)
            ms_mean = fitfun_poly5(mh_grid,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])       
            for i in range(mh_grid.size):
                num_ms += num_mhalo[i] * 1./(ms_std_grid[i]*np.sqrt(2.*np.pi)) * \
                    np.exp(-0.5 * ((ms_grid-ms_mean[i])/ms_std_grid[i])**2.)*dmh

            pylab.axis([6.,11.3,-5.,1.])
            ok=(ms_grid > np.min(ms)) & (ms_grid < np.max(ms))
            pylab.plot(ms_grid[ok],np.log10(num_ms[ok]),'k')

            fname='/Users/phopkins/Documents/work/plots/python/behroozi_sfh/stellar_mf/z'
            if (z_redshift==0.0): zkeys=['0.1']
            if (z_redshift==0.5): zkeys=['0.5']
            if (z_redshift==1.0): zkeys=['1.0']
            if (z_redshift==2.0): zkeys=['2.0']
            if (z_redshift==4.0): zkeys=['4.0']
            if (z_redshift==6.0): zkeys=['6.0']#,'7.0','8.0']
            for zkey in zkeys:
                v = np.loadtxt(fname+zkey+'.dat')
                x=v[:,0]; y=v[:,1]; dy=np.sqrt(v[:,2]**2.+v[:,3]**2.+0.05**2.)
                #pylab.plot(x,y,'mo')
                pylab.plot(x,y+dy,'m--')
                pylab.plot(x,y-dy,'m--')


        
        if (SHOW_MGAL_SCATTER):
            pylab.axis([8.5,13.5,0.,0.53])
            pylab.axis([8.5,12.,-0.02,0.55])
            pylab.xlabel(r'$\log(\ {\rm M}_{\rm halo}\ /\ {\rm M}_{\odot}\ )$',fontsize=charsi)
            pylab.ylabel(r'$\sigma[{\log({\rm M}_{\ast}\,|\,{\rm M}_{\rm halo})}]$')
            #pylab.axis([3.,11.2,0.,0.5])
        
            ok=(all_ms > 2.); ms=all_ms[ok]; mh=all_mh[ok]; 
            s=np.argsort(mh); mh=mh[s]; ms=ms[s];
            dmh=0.05; dms=1.*dmh; 
            dmh=0.25; dms=1.*dmh; 

            mh_grid = np.arange(mh[0],mh[-1],dmh); 
            #popt,pcov = optimize.curve_fit(fitfun_poly5,mh,ms,p0=[1.,1.,1.,1.,1.,1.]);
            #ms_expected = fitfun_poly5(mh,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])   
            #ms_grid = fitfun_poly5(mh_grid,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) 
            popt,pcov = optimize.curve_fit(fitfun_poly3,mh,ms,p0=[1.,1.,1.,1.]);
            ms_expected = fitfun_poly3(mh,popt[0],popt[1],popt[2],popt[3])   
            ms_grid = fitfun_poly3(mh_grid,popt[0],popt[1],popt[2],popt[3]) 
            ms_std = 0.*mh;
            d_i = 5
            for i in range(mh.size):
                imin=i-d_i; 
                if (imin<=0): imin=0
                imax=i+d_i; 
                if (imax>=mh.size): imax=mh.size
                d_ms_in=ms[imin:imax] - ms_expected[imin:imax]
                ms_std[i] = np.std(d_ms_in)/np.sqrt(2.)
            ms_std_grid = np.interp(mh_grid,mh,ms_std)
                
            #pylab.plot(mh,ms_std,colors_vec[i_z]+'o')
            pylab.plot(mh_grid,ms_std_grid,'ko:')
            ##pylab.plot(ms,ms_std,colors_vec[i_z]+'.')
            ##pylab.plot(ms_grid,ms_std_grid,colors_vec[i_z])
                
            s=np.argsort(ms); ms=ms[s]; mh=mh[s];
            ms_grid = np.arange(ms[0],ms[-1],dmh); 
            #popt,pcov = optimize.curve_fit(fitfun_poly5,ms,mh,p0=[1.,1.,1.,1.,1.,1.]);
            #mh_expected = fitfun_poly5(ms,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])       
            #mh_grid = fitfun_poly5(ms_grid,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) 
            popt,pcov = optimize.curve_fit(fitfun_poly3,ms,mh,p0=[1.,1.,1.,1.]);
            mh_expected = fitfun_poly3(ms,popt[0],popt[1],popt[2],popt[3])       
            mh_grid = fitfun_poly3(ms_grid,popt[0],popt[1],popt[2],popt[3]) 
            mh_std = 0.*ms;
            d_i = 20
            for i in range(ms.size):
                imin=i-d_i; 
                if (imin<=0): imin=0
                imax=i+d_i; 
                if (imax>=mh.size): imax=mh.size
                d_mh_in=mh[imin:imax] - mh_expected[imin:imax]
                mh_std[i] = np.std(d_mh_in-np.mean(d_mh_in)) /np.sqrt(2.)
            mh_std_grid = np.interp(ms_grid,ms,mh_std)
 
            #pylab.plot(mh,mh_std,colors_vec[i_z]+'o')
            #pylab.plot(mh_grid,mh_std_grid,colors_vec[i_z])
            ###pylab.plot(ms,mh_std,colors_vec[i_z]+'o')
            #pylab.plot(ms_grid,mh_std_grid,colors_vec[i_z])

        if((SHOW_MGAL_MHALO_ALTSIMS)|(SHOW_MGAL_MHALO_Z0)|(SHOW_MGAL_MHALO_NUMERICS)):
            pylab.subplots_adjust(left=(26./25.4)/8., 
                bottom = (17./25.4)/12., 
                right = 1 - (5./25.4)/8., 
                top = 1 - (2.5/25.4)/12.,
                hspace = 0.08)

        if(SHOW_MGAL_MHALO_ALLZ):
            pylab.subplots_adjust(left=0.09, bottom=0.082, right=1-0.05, top=1-0.07,
                hspace=0.12, wspace=0.10)

        if(SHOW_SSFR):
            pylab.subplots_adjust(left=0.09, bottom=0.082, right=1-0.05, top=1-0.07,
                hspace=0.26, wspace=0.12)
 
    if (SHOW_MGAL_MHALO_Z0):
        pylab.subplot(211)
        #pylab.legend(loc='lower right',fontsize=14,\
        #    frameon=False,labelspacing=0.2,numpoints=1)
        #pylab.legend([leg_09,leg_10,leg_11,leg_12v,leg_12q,leg_13],\
        #    ['m09','m10','m11','m12v','m12q','m13'],loc='lower right',\
        #    fontsize=14,frameon=False,labelspacing=0.3,numpoints=1,markerscale=0.75,\
        #    handletextpad=0.)
        pylab.legend([leg_09,leg_10,leg_11,leg_12v,leg_12i,leg_12q,leg_13],\
            ['m09','m10','m11','m12v','m12i','m12q','m13'],loc='lower right',\
            fontsize=14,frameon=False,labelspacing=0.3,numpoints=1,markerscale=0.75,\
            handletextpad=0.)
        pylab.savefig("mg_mh_z0.pdf",transparent=True,bbox_inches='tight',pad_inches=0)
    if (SHOW_MGAL_MHALO_ALTSIMS):
        pylab.savefig("mg_mh_altsims.pdf",transparent=True,bbox_inches='tight',pad_inches=0)
    if (SHOW_MGAL_MHALO_ALLZ):
        pylab.savefig("mg_mh_allz.pdf",transparent=True,bbox_inches='tight',pad_inches=0)
    if (SHOW_MGAL_MHALO_NUMERICS):
        pylab.subplot(211)
        pylab.legend(loc='upper left',fontsize=12,\
            frameon=False,labelspacing=0.2,numpoints=1)
        pylab.subplot(212)
        pylab.legend(loc='lower right',fontsize=14,\
            frameon=False,labelspacing=0.2,numpoints=1)
        pylab.savefig("mg_mh_fb.pdf",transparent=True,bbox_inches='tight',pad_inches=0)
    if (SHOW_MGAL_SCATTER):
        pylab.savefig("mg_mh_scatter.pdf",transparent=True,bbox_inches='tight',pad_inches=0)
    if (SHOW_SSFR):
        pylab.savefig("mg_mh_sfr.pdf",transparent=True,bbox_inches='tight',pad_inches=0)






def fitfun_poly3( x, p0, p1, p2, p3):
    return p0+p1*x+p2*x*x+p3*x*x*x

def fitfun_poly5( x, p0, p1, p2, p3, p4, p5):
    return p0+p1*x+p2*x*x+p3*x*x*x+p4*x*x*x*x+p5*x*x*x*x*x




def quicklook(sdir, snum, cen, rcut):
    pylab.close('all')
    pylab.axis([-rcut,rcut,-rcut,rcut]);
    for ptype in [3,2,1,0,4]:
        P=gadget.readsnap(sdir,snum,ptype,cosmological=1,skip_bh=1);
        if (P['k']==1): 
            if (P['m'].size>1):
                p=P['p']; x=p[:,0]-cen[0]; y=p[:,1]-cen[1]; z=p[:,2]-cen[2];
                r=np.sqrt(x*x+y*y+z*z);
                ok=(r<rcut)
                colors=['k','r','b','g','m','c','y']
                colors=['g','k','m','b','r','c','y']
                pylab.plot(x[ok],y[ok],colors[ptype]+',')
    


def ics_scan():
    sdir='/Users/phopkins/Documents/work/plots/zooms/m14_tst'
    snum=7
    infiname = halo_finder.get_halo_id_filename(sdir,snum,full_halo_results=1,four_char=0)
    infi=h5py.File(infiname,'r')
    
    mh=np.log10(infi['halo_m'])+10.
    
    print(mh[0:50])
    
    pos=infi['halo_pos']
    
    ok = (mh >= 12.)
    pylab.scatter(pos[ok,0],pos[ok,1],marker='.',color='k')

    ok = (mh >= 12.5)
    pylab.scatter(pos[ok,0],pos[ok,1],marker='.',color='y')

    ok = (mh >= 13.)
    pylab.scatter(pos[ok,0],pos[ok,1],marker='.',color='b')

    ok = (mh >= 13.5)
    pylab.scatter(pos[ok,0],pos[ok,1],marker=',',color='g')

    ok = (mh >= 14.)
    pylab.scatter(pos[ok,0],pos[ok,1],marker=',',color='m')

    print(infi['halo_pos'][ok,:])
    print(np.log10(infi['halo_m'][ok]*0.7)+10.)
    print(infi['halo_r_med'][ok])
    print(infi['halo_r_max'][ok])
    print(infi['halo_r_rms'][ok])

    ok = (mh >= 14.5)
    pylab.scatter(pos[ok,0],pos[ok,1],marker=',',color='r')

    
    infi.close()



def dm_profile(sdir, snum, cen):
	pos=np.zeros([1,3],dtype=float)
	m=np.zeros([1],dtype=float)
	for ptype in [1,2,3,5]:
		P=gadget.readsnap(sdir,snum,ptype,cosmological=1,skip_bh=1);
		if (P['k']==1): 
			if (ptype<5) or (P['m'].shape[0]>200):
				pos=np.concatenate((pos,P['p'])); 
				m=np.concatenate((m,P['m']));
	r=np.sqrt((pos[:,0]-cen[0])**2 + (pos[:,1]-cen[1])**2 + (pos[:,2]-cen[2])**2);
	s=np.argsort(r); s=s[1:s.shape[0]]; r=r[s]; m=np.cumsum(m[s]);
	ln_r=np.log(r); ln_m=np.log(m);

	dr=0.1;
	rg=np.arange(ln_r[0],ln_r[-1],dr);
	f_mg=interpolate.interp1d(ln_r,ln_m,kind='linear'); mg=f_mg(rg);
	dlnm_dlnr = np.diff(mg)/np.diff(rg);
	dlnm_dlnr = np.concatenate(([dlnm_dlnr[0]],dlnm_dlnr)) + 1.0e-10;

	log_rho = 10.+(mg - 3.*rg - np.log(4.*math.pi) + np.log(dlnm_dlnr))/np.log(10.);#
	log_r = rg/np.log(10.)
	
	pylab.axis([10.**(-1.4),10.**(2.3),0.7,8.7]); pylab.xscale('log');
	pylab.xlabel(r'$r\ \ [{\rm kpc}]$')
	pylab.ylabel(r'$\rho_{\rm DM}\ \  [{\rm M_{\odot}\,kpc^{-3}}]$')
	pylab.xticks([0.1,1.,10.,100.],['0.1','1','10','100'])
	#pylab.scatter(10.**log_r,log_rho,marker='o',facecolors='none',edgecolors='k')
	pylab.plot(10.**log_r,log_rho,linewidth=1.0,color='k')
	
	ok = (-1. < log_r) & (log_r < 2.)
	popt,pcov = optimize.curve_fit(fitfun.nfw,log_r[ok],log_rho[ok],p0=[7.,0.5]);
	pylab.plot(10.**log_r,fitfun.nfw(log_r,popt[0],popt[1]),'r:',linewidth=2.0);
	print('dm halo fit parameters == '); print(popt)

	pylab.savefig("dm_profile.ps"); pylab.close('all')


def vc_profile(sdir, snum, cen):
	pos=np.zeros([1,3],dtype=float)
	m=np.zeros([1],dtype=float)
	for ptype in [0,1,2,3,4,5]:
		P=gadget.readsnap(sdir,snum,ptype,cosmological=1,skip_bh=1);
		if (P['k']==1): 
			if (ptype<5) or (P['m'].shape[0]>200):
				pos=np.concatenate((pos,P['p'])); 
				m=np.concatenate((m,P['m']));
	r=np.sqrt((pos[:,0]-cen[0])**2 + (pos[:,1]-cen[1])**2 + (pos[:,2]-cen[2])**2);
	s=np.argsort(r); s=s[1:s.shape[0]]; r=r[s]; m=np.cumsum(m[s]); ln_r=np.log(r);

	dr=0.1;
	vc = np.sqrt(6.67e-8 * (m*2.0e43) / (3.086e21*r)) / 1.0e5; # km/s
	rg=np.arange(ln_r[0],ln_r[-1],dr);
	f_vc=interpolate.interp1d(ln_r,np.log(vc),kind='linear'); vc_g=np.exp(f_vc(rg));
	
	pylab.axis([10.**(-1.4),10.**(2.3), 0.,350.])
	pylab.xticks([0.1,1.,10.,100.],['0.1','1','10','100']); pylab.xscale('log')
	pylab.axis([0.,40., 0.,50.]); 
	pylab.axis([0.,40., 0.,200.]); 
	pylab.axis([0.,60., 0.,350.]); 
	pylab.xscale('linear')
	pylab.xlabel(r'$r\ \ [{\rm kpc}]$')
	pylab.ylabel(r'$V_{c}\ \  [{\rm km\,s^{-1}}]$')
	pylab.plot(np.exp(rg),vc_g,linewidth=1.0,color='k')

	pylab.savefig("vc_profile.ps"); pylab.close('all')
	
		

def sfr_time(sdir, snum, cen, rcut=10.):
	dz=0.005
	P=gadget.readsnap(sdir,snum,4,cosmological=1,skip_bh=1);
	a_form=P['age']; m=P['m'];
	s=np.argsort(a_form); a_form=a_form[s]; m=np.cumsum(m[s]);
	z_grid=10.**np.arange(0.,np.log10(21.),dz)-1.; a_grid=1./(1.+z_grid);
	t_grid=cosmo.lookback_time(z_grid).value;
	ok=(a_grid >= a_form[0]) & (a_grid <= a_form[-1]);
	z_grid=z_grid[ok]; a_grid=a_grid[ok]; t_grid=t_grid[ok];
	
	f_mg=interpolate.interp1d(a_form,m,kind='linear'); mg=f_mg(a_grid);
	dm_dt=-np.diff(mg)/np.diff(t_grid);
	t_mid=t_grid[0:t_grid.shape[0]-1]+np.diff(t_grid); # midpoints where the derivatives apply
	z_mid=z_grid[0:z_grid.shape[0]-1]+np.diff(z_grid);
	dm_dt *= (1.e10/1.e9) / 0.7; # units and mass loss
	print(np.max(dm_dt))
	pylab.close('all'); 
	#pylab.axis([0.,1., 0., 0.17]); pylab.xscale('linear'); pylab.yscale('linear')
	pylab.axis([0.,12., 0., np.max(dm_dt)*1.1]); pylab.xscale('linear'); pylab.yscale('linear')
	#pylab.axis([0.,12., 0.01, np.max(dm_dt)*1.1]); pylab.xscale('linear'); pylab.yscale('log')
	pylab.xlabel('Redshift z')
	pylab.ylabel(r'${\rm SFR}\ \ \dot{M}_{\ast}\ \  [{\rm M_{\odot}\,yr^{-1}}]$')
	pylab.plot(z_mid,dm_dt,linewidth=1.0,color='k')
	pylab.savefig("sfr_hist.ps");
	pylab.close('all');



def sb_profile(sdir, snum, cen):
    P=gadget.readsnap(sdir,snum,4,cosmological=1,skip_bh=1);
    P_head=gadget.readsnap(sdir,snum,4,cosmological=1,skip_bh=1,header_only=1);
    m=P['m']*1.0e10; zm=P['z']; pos=P['p']; z_eff=zm[:,0]/0.02;
    x=pos[:,0]-cen[0]; y=pos[:,1]-cen[1]; z=pos[:,2]-cen[2];
    r=np.sqrt(x*x+y*y+z*z)
    m1=math.fsum(m[r<1.])
    m3=math.fsum(m[r<3.])
    m5=math.fsum(m[r<5.])
    m10=math.fsum(m[r<10.])
    m20=math.fsum(m[r<20.])
    print('mass inside r = 1, 3, 5, 10, 20:',np.log10(m1),np.log10(m3),np.log10(m5),np.log10(m10),np.log10(m20))
    #z_min=0.01; z_max=3.; z_eff[z_eff < z_min] = z_min; z_eff[z_eff > z_max] = z_max;

    dlnr=0.1*np.log(10.);
    BAND=2;
    r=np.sqrt(x*x+y*y);

    age=gadget.get_stellar_ages(P,P_head,cosmological=1);
    lum = m * sps.colors_table(age, z_eff, BAND_ID=BAND, CRUDE=0);
    s=np.argsort(r); ln_r=np.log(r[s]); m=np.cumsum(m[s]); lum=np.cumsum(lum[s]); 
    ln_rg=np.arange(ln_r[0],np.log(100.),dlnr);
    
    # mass surface density
    lnx=np.interp(ln_rg,ln_r,np.log(m)); 
    xmid=lnx[0:lnx.shape[0]-1]+np.diff(lnx); rmid=ln_rg[0:ln_rg.shape[0]-1]+np.diff(ln_rg);
    sigma_m=np.exp(xmid)/(2.*math.pi*np.exp(2.*rmid)) * np.diff(lnx)/np.diff(ln_rg);

    # luminosity surface density
    lnx=np.interp(ln_rg,ln_r,np.log(lum)); 
    xmid=lnx[0:lnx.shape[0]-1]+np.diff(lnx); rmid=ln_rg[0:ln_rg.shape[0]-1]+np.diff(ln_rg);
    sigma_l=np.exp(xmid)/(2.*math.pi*np.exp(2.*rmid)) * np.diff(lnx)/np.diff(ln_rg);

    pylab.close('all'); 
    pylab.axis([0.01, 100., 4., 11.]); pylab.xscale('log')
    pylab.xlabel(r'$R\ \ [{\rm kpc}]$')
    pylab.ylabel(r'$\log(\Sigma_{\ast})\ \  [{\rm M_{\odot}\,kpc^{-2}}]$')
    rr=np.exp(rmid)
    rr=rr**0.25; pylab.xscale('linear'); pylab.axis([0.,3.2,4.,11.])
    pylab.plot(rr,np.log10(sigma_m),'k',linewidth=1.0)
    pylab.plot(rr,np.log10(sigma_l),'o--',linewidth=1.0)

    x=rmid/np.log(10.); ym=np.log10(sigma_m); yl=np.log10(sigma_l);
    #ok=(x > -2.) & (x < 1.0) & (ym > 2.0);
    ok=(x > -2.) & (x < 1.6) & (ym > 2.0);
    #ok=(x > -2.) & (x < 0.5) & (ym > 2.0);

    popt,pcov = optimize.curve_fit(fitfun.sersic,x[ok],ym[ok],p0=[8.,0.,1.]);
    pylab.plot(rr,fitfun.sersic(x,popt[0],popt[1],popt[2]),'r:',linewidth=2.0);
    print('(mass) sersic fit == '); print(popt); popt_mass=popt;
    popt,pcov = optimize.curve_fit(fitfun.sersic_disk,x[ok],ym[ok],p0=[8.,1.,1., 9.,-0.5,1.]);
    yfit = fitfun.sersic_disk(x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
    pylab.plot(rr,yfit,'m-.',linewidth=2.0);
    print('(mass) sersic-disk fit == '); print(popt); popt_sd_mass=popt;
    ## integrate this model
    #dm=yfit*2.*np.pi*rr*rr*dlnr
    

    popt,pcov = optimize.curve_fit(fitfun.sersic,x[ok],yl[ok],p0=popt_mass);
    pylab.plot(rr,fitfun.sersic(x,popt[0],popt[1],popt[2]),'g:',linewidth=2.0);
    print('(lum) sersic fit == '); print(popt)
    popt,pcov = optimize.curve_fit(fitfun.sersic_disk,x[ok],yl[ok],p0=popt_sd_mass);
    pylab.plot(rr,fitfun.sersic_disk(x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),'c-.',linewidth=2.0);
    print('(lum) sersic-disk fit == '); print(popt)
    
    pylab.savefig("sb_profile.ps");
    pylab.close('all');
    
    
def ang_mom(sdir, snum, cen):
    # first get the Vc profile :: 
    pos=np.zeros([1,3],dtype=float)
    m=np.zeros([1],dtype=float)
    for ptype in [0,1,2,3,4,5]:
        P=gadget.readsnap(sdir,snum,ptype,cosmological=1,skip_bh=1);
        if (P['k']==1): 
            if (ptype<5) or (P['m'].shape[0]>200):
                pos=np.concatenate((pos,P['p'])); 
                m=np.concatenate((m,P['m']));
    r=np.sqrt((pos[:,0]-cen[0])**2 + (pos[:,1]-cen[1])**2 + (pos[:,2]-cen[2])**2);
    s=np.argsort(r); s=s[1:s.shape[0]]; r=r[s]; m=np.cumsum(m[s]); ln_r=np.log(r);
    
    dr=0.1;
    vc = np.sqrt(6.67e-8 * (m*2.0e43) / (3.086e21*r)) / 1.0e5; # km/s
    rg=np.arange(ln_r[0],ln_r[-1],dr);
    f_vc=interpolate.interp1d(ln_r,np.log(vc),kind='linear'); vc_g=np.exp(f_vc(rg));
    rg_vc=np.exp(rg); vc_rg=vc_g;
    
    # now load the stellar data
    P=gadget.readsnap(sdir,snum,4,cosmological=1,skip_bh=1);
    P_head=gadget.readsnap(sdir,snum,4,cosmological=1,skip_bh=1,header_only=1)
    pos=P['p']; xs=pos[:,0]-cen[0]; ys=pos[:,1]-cen[1]; zs=pos[:,2]-cen[2];
    vel=P['v']; vx=vel[:,0]; vy=vel[:,1]; vz=vel[:,2];
    ms=P['m']; rs=np.sqrt(xs*xs+ys*ys+zs*zs);
    zm=P['z']; z_eff=zm[:,0]/0.02;
    age=gadget.get_stellar_ages(P,P_head,cosmological=1);
    lum = ms * sps.colors_table(age,z_eff,BAND_ID=2,CRUDE=1);
    
    rmin=0.; rmax=3.;
    rmin=3.; rmax=10.;
    ok=(rs > rmin) & (rs < rmax) & (age < 3.0);
    ok=(rs > rmin) & (rs < rmax) & (age < 8.0);
    m=ms[ok]; lum=lum[ok]; x=xs[ok]; y=ys[ok]; z=zs[ok]; vx=vx[ok]; vy=vy[ok]; vz=vz[ok]
    wt=m/math.fsum(m); com_vx=math.fsum(wt*vx); com_vy=math.fsum(wt*vy); com_vz=math.fsum(wt*vz);
    
    vx=vx-com_vx; vy=vy-com_vy; vz=vz-com_vz; 
    r=np.sqrt(x*x+y*y+z*z);
    jx=(vy*z-vz*y); jy=-(vx*z-vz*x); jz=(vx*y-vy*x)
    vc_r = np.interp(r,rg_vc,vc_rg); jc_r = vc_r * r;

    jxm=math.fsum(jx); jym=math.fsum(jy); jzm=math.fsum(jz);
    jmm=np.sqrt(jxm*jxm+jym*jym+jzm*jzm);
    jxm/=jmm; jym/=jmm; jzm/=jmm;
    j_mag=np.sqrt(jx*jx+jy*jy+jz*jz);
    j_pro=jx*jxm+jy*jym+jz*jzm;

    jp_jc = j_pro / jc_r;

    xg=np.arange(-10.,10.,0.05); xm=xg[0:xg.shape[0]-1]+np.diff(xg);
    s=np.argsort(jp_jc); js=jp_jc[s]; m_s=np.cumsum(m[s]); l_s=np.cumsum(lum[s]);
    m_s=m_s/max(m_s); l_s=l_s/max(l_s);
    u=m_s; y=np.interp(xg,js,u); dm_dj=np.diff(y)/np.diff(xg);
    u=l_s; y=np.interp(xg,js,u); dl_dj=np.diff(y)/np.diff(xg);
    
    pylab.close('all'); #pylab.savefig("ang_mom.ps");
    pylab.axis([-2., 2., 0., 1.6]); pylab.xscale('linear')
    pylab.xlabel(r'$j_{z}\ /\ j_{c}$')
    pylab.ylabel(r'${\rm Fraction}$')
    pylab.plot(xm,dm_dj,'r:',linewidth=1.0)
    pylab.plot(xm,dl_dj,'k',linewidth=2.0)
    

def calculate_zoom_center(sdir,snum,cen=[0.,0.,0.],clip_size=2.e10):
    rgrid=np.array([1.0e10,1000.,700.,500.,300.,200.,100.,70.,50.,30.,20.,10.,5.,2.5,1.,0.5]);
    rgrid=rgrid[rgrid <= clip_size];
    Ps=gadget.readsnap(sdir,snum,4,cosmological=1);
    n_new=Ps['m'].shape[0];
    if (n_new > 1):
        pos=Ps['p']; x0s=pos[:,0]; y0s=pos[:,1]; z0s=pos[:,2];
    Pg=gadget.readsnap(sdir,snum,0,cosmological=1);
    rho=np.array(Pg['rho'])*407.5;
    if (rho.shape[0] > 0):
        pos=Pg['p']; x0g=pos[:,0]; y0g=pos[:,1]; z0g=pos[:,2];
    rho_cut=1.0e-5;
    cen=np.array(cen);

    for i_rcut in range(len(rgrid)):
        for j_looper in range(5):
            if (n_new > 1000):
                x=x0s; y=y0s; z=z0s;
            else:
                ok=(rho > rho_cut);
                x=x0g[ok]; y=y0g[ok]; z=z0g[ok];
            x=x-cen[0]; y=y-cen[1]; z=z-cen[2];
            r = np.sqrt(x*x + y*y + z*z);
            ok = (r < rgrid[i_rcut]);
            if (len(r[ok]) > 1000):
                x=x[ok]; y=y[ok]; z=z[ok];
                if (i_rcut <= len(rgrid)-5):
                    cen+=np.array([np.median(x),np.median(y),np.median(z)]);
                else:
                    cen+=np.array([np.mean(x),np.mean(y),np.mean(z)]);
            else:
                if (len(r[ok]) > 200):
                    x=x[ok]; y=y[ok]; z=z[ok];
                    cen+=np.array([np.mean(x),np.mean(y),np.mean(z)]);
                    
    return cen;



def sfh_z0_burstiness():
    pylab.close('all'); 
    sdir0='/Users/phopkins/Documents/work/plots/zooms/'
    bdir0='/Users/phopkins/Documents/work/plots/python/behroozi_sfh/sfh_individual/'

    plot.figure(1,figsize=(8.,6.))
    charsi=20.
    matplotlib.rcParams.update({'font.size':18})

    sdirs=['m09_dw' ,'m10_dw','m11_lr','hires_bh','m12_mr','m13_ics']
    snums=[ 339,      440   ,   440,    441    ,    441    ,   66  ]
    bspecs=['-1',      '-1' ,   '11',   '12',       '12'    ,   '13']
    ymins=[ 1.e-4   ,  0.00005, 0.03,   0.1   ,      0.3   ,   0.3]
    ymaxs=[ 0.02,        0.02,    3.,    9.   ,    40.     ,   40.]
    rcuts=[ 15.,          15.,    3.,      20.  ,    5.     ,    50.]        

    sdirs=['m09_dw','m10_dw','m11_lr','hires_bh','m12_mr','m13_mr']
    snums=[ 339,      440   ,   440,    441    ,    441    ,   437  ]
    bspecs=['-1',      '-1' ,   '11',   '12',       '12'    ,   '13']
    ymins=[ 1.e-4   ,  0.00005, 0.03,   0.1   ,      0.3   ,   0.3]
    ymaxs=[ 0.02,        0.02,    3.,    9.   ,    40.     ,   100.]
    rcuts=[ 15.,          15.,    3.,      20.  ,    5.     ,    50.]        
    slbls=['m09',    'm10',    'm11',   'm12v',   'm12q',  'm13']

    sdirs=['m09_dw','m10_hr_Aug8','m11_hr_Aug8','B1_Aug27_egylim','m12q_hr_Sep10_nolim','m13_mr_Aug8']
    snums=[ 339,      440   ,   440,    440    ,    440    ,   427  ]
    x000=  [2.875,3.2005, 42.1925,    7.0343, 39.668,     40.3695]
    y000=  [2.2626, 2.9754, 45.8462,    1.3945, 48.5812,    46.0468]
    z000=  [2.1408,  3.2835, 40.18,      9.417,  44.851,     34.234 ]
    bspecs=['-1',      '-1' ,   '11',   '-1',       '12'    ,   '13']
    ymins=0.1*np.array([ 1.e-4   ,  0.00005, 0.03,   0.1   ,      0.3   ,   0.3])
    ymaxs=[ 0.02,        0.02,    3.,    9.   ,    40.     ,   100.]
    rcuts=np.array([ 100.,          50.,    150.,      200.  ,    400.     ,    500.])  * 10.
    #rcuts=1.e10*np.array(rcuts)
    slbls=['m09',    'm10',    'm11',   'm12v',   'm12q',  'm13']

    marker_syms=['*','s','D','o','p','x','D']
    colors_vec=[ 'c','b','g','k','r','m','y']

    sdirs=['m09_Jan2014','m10_hr_Jan2014','m11_hr_Jan2014','B1_11_Jan2014','m12v_3_Jan2014',
            'm12qq_Jan2014','m13_mr_Jan2014']#,'m14_Jan2014']
    snums=[ 440,      440   ,   440,    440    ,    440    ,   440, 440, 440  ]
    x000=  [4.706,  3.207635,42.051,6.97123,41.92507,39.48788,43.0227,44.6275]
    y000=  [3.42459,2.976621,45.901, 1.44426,44.67706,47.9995, 49.3577,46.4436]
    z000=  [3.17808,3.286555,39.757,9.42507,46.069,  44.80903,36.1775,37.0167]
    bspecs=['-1',      '-1' ,   '11',   '12',       '12' ,'12'    ,   '13','14']
    ymins=0.1*np.array([ 1.e-5   ,  0.00005, 0.03,   0.3   ,      0.3   ,  0.3, 0.3, 3.])
    ymaxs=[ 0.002,        0.02,    3.,    40.   ,    40.    , 40. ,   100., 1000.]
    rcuts=[ 20.,          20.,    100.,      40.  ,    40.     ,  500.,  500.,   200.]   
    #rcuts=1.e10*np.array(rcuts)
    slbls=['m09',    'm10',    'm11',   'm12v',   'm12i',  'm12q', 'm13', 'm14']

    marker_syms_0=np.array(['o','s','p','^','*','x','D','+','h','1','v','2','<','3','>','4','H'])
    colors_vec_0=np.array(['black','blue','red','green','deepskyblue','darkviolet','orange','magenta','gold','sienna','pink','forestgreen','darkgray'])
    i_order=[4,1,3,0,6,2,5,7]
    marker_syms=marker_syms_0[i_order]
    colors_vec=colors_vec_0[i_order]

    #sdirs=['m11_lr']; snums=[440]; bspecs=['11']; ymins=[0.03]; ymaxs=[3.]; rcuts=[3.];

    dt_grid = np.arange(6.,10.,0.3)
    dt_units = (10.**(dt_grid)) / 1.e9
    dt_0 = dt_units[0]
    dt_rescale = np.around(dt_units/dt_0/2.)
    dt_rescale = 2*dt_rescale.astype(int)+1

    for i,sdir,snum,bspec,ymin,ymax,rcut in \
      zip(range(np.array(sdirs).size),sdirs,snums,bspecs,ymins,ymaxs,rcuts):
        if (i>-1):
            print(sdir)
            z,sfr = sfr_time_snap( sdir0+sdir, snum, rcut=rcut, dt=dt_0 ,\
                             cen=[x000[i],y000[i],z000[i]])
            mcorr_mloss = 0.7 + 0.3/(1.+z)
            sfr /= mcorr_mloss
            xx = 1.+z
            zmax = 12.
            #pylab.axis([1.+0.,1.+zmax,ymin,ymax])
            #pylab.xscale('log'); pylab.yscale('log')
            #ztick=1.+np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
            #pylab.xticks(ztick,['0','1','2','3','4','5','6','7','8','','10','',''])
            #pylab.xlabel(r'Redshift $z$')
            #pylab.ylabel(r'${\rm SFR}\ \ \dot{M}_{\ast}\ \  [{\rm M_{\odot}\,yr^{-1}}]$')
            
            #pylab.axis([6.,10.,0.,1.])
            pylab.axis([6.,10.,-0.025,1.025])
            pylab.xlabel(r'$\log(\Delta t_{\rm avg}\ /\ {\rm yr})$')
            pylab.ylabel(r'Scatter $\sigma_{\rm SFR}$  [dex]')
            pylab.xscale('linear'); pylab.yscale('linear')
            
            s=np.argsort(xx); xx=xx[s]; sfr=sfr[s];

            #sfr_smoo_70 = sfr
            #sfr_smoo_75 = util.smooth(sfr_smoo_70, window_len=3, window='flat')
            #sfr_smoo_80 = util.smooth(sfr_smoo_70, window_len=11, window='flat')
            #sfr_smoo_85 = util.smooth(sfr_smoo_70, window_len=33, window='flat')
            #sfr_smoo_90 = util.smooth(sfr_smoo_70, window_len=111, window='flat')
            #sfr_smoo_95 = util.smooth(sfr_smoo_70, window_len=333, window='flat')
            #sfr_smoo_100= util.smooth(sfr_smoo_70, window_len=1111, window='flat')

            #pylab.plot(xx, sfr_smoo_70, 'b:')
            #pylab.plot(xx, sfr_smoo_80, 'k--')
            #pylab.plot(xx, sfr_smoo_90, 'r-')
            #pylab.plot(xx, sfr_smoo_95, 'g-', linewidth=2)
            
            winsi=np.copy(dt_rescale)
            ii = np.arange(0.,xx.size,1.)
            scatter=0.*winsi
            for windowsize,j in zip(winsi,range(winsi.size)):
                if (sfr.size > windowsize):
                    if (windowsize<=1):
                        sfr_smoo = 1.*sfr;
                    else:
                        sfr_smoo = util.smooth(sfr, window_len=windowsize, window='flat');
                    ok = (np.isnan(sfr_smoo)==False) & \
                         (sfr_smoo > 0.01*ymin) & \
                         (sfr_smoo < 100.*ymax) & \
                         (xx > 0.) & (xx < 10.) & \
                         (ii > np.float(windowsize)/2.) & \
                         (ii < (np.float(xx.size)-np.float(windowsize)/8.)) 
                    s_ok = sfr_smoo[ok]
                    if (s_ok.size > 10):
                        x_ok = xx[ok]
                        #pylab.plot(x_ok, s_ok)
                        xg = np.arange(np.log(np.min(x_ok))+0.025,np.log(np.max(x_ok))-0.025,\
                                np.log(np.max(x_ok)/np.min(x_ok))/200.)
                        yg = np.interp(xg, np.log(x_ok), np.log(s_ok))
                        ndeg = 8
                        ndeg = 12
                        Pfit = np.polyfit(xg,yg,ndeg)
                        y_fit = 0.*x_ok
                        for k in range(ndeg+1):
                            y_fit += Pfit[Pfit.size-k-1]*(np.log(x_ok)**np.float(k))
                        #pylab.plot(np.exp(xg), np.exp(yg), ':')
                        #pylab.plot(x_ok, np.exp(y_fit), '--')
                        delta = np.log(s_ok) - y_fit
                        scatter[j] = np.std(delta)
            scatter_dex = scatter / np.log(10.)
            print(scatter_dex)
            
            markersizer=11.
            if(marker_syms[i]=='*'): 
                markersizer=18.
            pylab.plot(dt_grid,scatter_dex,color=colors_vec[i],linestyle=':')
            pylab.plot(dt_grid,scatter_dex,marker=marker_syms[i],
                markerfacecolor='w',markeredgecolor=colors_vec[i],
                color=colors_vec[i],markersize=markersizer,linestyle='',
                markeredgewidth=2,label=slbls[i])

    pylab.legend(loc='upper right',\
        fontsize=14,frameon=False,labelspacing=0.3,numpoints=1,markerscale=0.75,\
        handletextpad=0.)
    pylab.savefig("zoom_sfh_bursty.pdf",transparent=True,bbox_inches='tight',pad_inches=0)

# m09_dw : cen=[2.875,2.2626,2.1408]
# m09_Aug8 : cen=[4.70515,3.4245,3.178]
# m10_dw : cen=[3.3066,3.0296,3.2355]
# m10_hr_Aug8 : cen=[3.2005,2.9754,3.2835]
# m11_hr_Aug8 : cen=[42.1925,45.8462,40.18]
# m11_lr : cen=[6.599,5.647,3.5025]
# hires_bh : cen=[7.027,1.424,9.416]
# B1_Aug27_egylim : cen=[7.0343,1.3945,9.417]
# m12q_mr_Aug8 : cen=[38.288,47.853,45.054]
# m12q_hr_Sep10_nolim : cen=[39.668,48.5812,44.851]
# m13_mr : cen=[41.9963,46.8648,37.1041]
# m13_mr_Aug8 : cen=[40.3695,46.0468,34.234]


def sfh_z0_comparison():
    pylab.close('all'); 
    sdir0='/Users/phopkins/Documents/work/plots/zooms/'
    bdir0='/Users/phopkins/Documents/work/plots/python/behroozi_sfh/sfh_individual/'

    #plot.figure(1,figsize=(24.,6.))
    plot.figure(1,figsize=(24.,12.))
    #plot.figure(1,figsize=(24.,12.))
    charsi=20.
    matplotlib.rcParams.update({'font.size':18})

    sdirs=['m09_dw','m10_dw','m11_lr','hires_bh','m12_mr','m13_ics']
    snums=[ 339,      440   ,   440,    441    ,    441    ,   66  ]
    bspecs=['-1',      '-1' ,   '11',   '12',       '12'    ,   '13']
    ymins=[ 1.e-4   ,  0.00005, 0.03,   0.1   ,      0.3   ,   0.3]
    ymaxs=[ 0.02,        0.02,    3.,    9.   ,    40.     ,   40.]
    rcuts=[ 15.,          15.,    3.,      20.  ,    5.     ,    50.]        
    slbls=['m09',    'm10',    'm11',   'm12v',   'm12q',  'm13']

    sdirs=['m09_dw','m10_dw','m11_lr','hires_bh','m12_mr','m13_mr']
    snums=[ 339,      440   ,   440,    441    ,    441    ,   437  ]
    x000= [2.875,   3.3066, 6.599,  7.027,  38.288, 41.9963]
    y000= [2.2626,  3.0296, 5.647,  1.424,  47.853, 46.8648]
    z000= [2.1408,  3.2355, 3.5025, 9.416,  45.054, 37.1041]
    bspecs=['-1',      '-1' ,   '11',   '12',       '12'    ,   '13']
    ymins=[ 1.e-4   ,  0.00005, 0.03,   0.1   ,      0.3   ,   0.3]
    ymaxs=[ 0.02,        0.02,    3.,    9.   ,    40.     ,   100.]
    rcuts=[ 15.,          15.,    3.,      20.  ,    20.     ,    50.]        
    slbls=['m09',    'm10',    'm11',   'm12v',   'm12q',  'm13']

    sdirs=['m09_Aug8','m10_hr_Aug8','m11_hr_Aug8','B1_Aug27_egylim','m12q_hr_Sep10_nolim','m13_mr_Aug8']
    snums=[ 441,      440   ,   440,    440    ,    440    ,   427  ]
    x000=  [4.70515,3.2005, 42.1925,    7.0343, 39.668,     40.3695]
    y000=  [3.4245, 2.9754, 45.8462,    1.3945, 48.5812,    46.0468]
    z000=  [3.178,  3.2835, 40.18,      9.417,  44.851,     34.234 ]
    bspecs=['-1',      '-1' ,   '11',   '12',       '12'    ,   '13']
    ymins=[ 1.e-4   ,  0.00005, 0.03,   0.1   ,      0.3   ,   0.3]
    ymaxs=[ 0.02,        0.02,    3.,    9.   ,    40.     ,   100.]
    rcuts=[ 20.,          5.,    15.,      20.  ,    20.     ,    15.]        
    slbls=['m09',    'm10',    'm11',   'm12v',   'm12q',  'm13']

    sdirs=['m09_dw','m10_hr_Aug8','m11_hr_Aug8','B1_Aug27_egylim','m12q_hr_Sep10_nolim','m13_mr_Aug8']
    snums=[ 339,      440   ,   440,    440    ,    440    ,   427  ]
    x000=  [2.875,3.2005, 42.1925,    7.0343, 39.668,     40.3695]
    y000=  [2.2626, 2.9754, 45.8462,    1.3945, 48.5812,    46.0468]
    z000=  [2.1408,  3.2835, 40.18,      9.417,  44.851,     34.234 ]
    bspecs=['-1',      '-1' ,   '11',   '-1',       '12'    ,   '13']
    ymins=0.1*np.array([ 1.e-4   ,  0.00005, 0.03,   0.1   ,      0.3   ,   0.3])
    ymaxs=[ 0.02,        0.02,    3.,    9.   ,    40.     ,   100.]
    ymins=0.1*np.array([ 1.e-4   ,  0.00005, 0.011,   0.05   ,      0.1   ,   0.3])
    ymaxs=[ 0.02,        0.02,    3.,    9.   ,    40.     ,   300.]
    rcuts=[ 20.,          3.,    20.,      20.  ,    40.     ,    100.]   
    #rcuts=1.e10*np.array(rcuts)
    slbls=['m09',    'm10',    'm11',   'm12v',   'm12q',  'm13']

    sdirs=['m09_Jan2014','m10_lr_Jan2014','m11_lr_Jan2014','B1_11_Jan2014','m12v_3_Jan2014',
            'm12qq_Jan2014','m13_mr_Jan2014']#,'m14_Jan2014']
    snums=[ 440,      440   ,   440,    440    ,    440    ,   440, 440, 440  ]
    x000=  [4.706,  3.3083,42.222,6.97123,41.92507,39.48788,43.0227,44.6275]
    y000=  [3.42459,3.0315,45.83, 1.44426,44.67706,47.9995, 49.3577,46.4436]
    z000=  [3.17808,3.2387,40.165,9.42507,46.069,  44.80903,36.1775,37.0167]
    bspecs=['-1',      '-1' ,   '11',   '12',       '12' ,'12'    ,   '13','14']
    ymins=0.1*np.array([ 1.e-5   ,  0.00005, 0.03,   0.3   ,      0.3   ,  0.3, 0.3, 3.])
    ymaxs=[ 0.002,        0.02,    3.,    40.   ,    40.    , 40. ,   100., 1000.]
    rcuts=[ 20.,          3.,    20.,      40.  ,    40.     ,  100.,  300.,   200.]   
    #rcuts=1.e10*np.array(rcuts)
    slbls=['m09',    'm10',    'm11',   'm12v',   'm12i',  'm12q', 'm13', 'm14']

    if(2==0):
        sdirs=['m09_Jan2014','m10_hr_Jan2014','m11_hr_Jan2014','B1_11_Jan2014','m12v_3_Jan2014',
                'm12qq_Jan2014','m13_mr_Jan2014']#,'m14_Jan2014']
        snums=[ 440,      440   ,   440,    440    ,    440    ,   440, 440, 440  ]
        x000=  [4.706,  3.207635,42.051,6.97123,41.92507,39.48788,43.0227,44.6275]
        y000=  [3.42459,2.976621,45.901, 1.44426,44.67706,47.9995, 49.3577,46.4436]
        z000=  [3.17808,3.286555,39.757,9.42507,46.069,  44.80903,36.1775,37.0167]
        bspecs=['-1',      '-1' ,   '11',   '12',       '12' ,'12'    ,   '13','14']
        ymins=0.1*np.array([ 1.e-5   ,  0.00005, 0.03,   0.3   ,      0.3   ,  0.3, 0.3, 3.])
        ymaxs=[ 0.002,        0.02,    3.,    40.   ,    40.    , 40. ,   100., 1000.]
        rcuts=[ 20.,          20.,    100.,      40.  ,    40.     ,  500.,  500.,   200.]   
        #rcuts=1.e10*np.array(rcuts)
        slbls=['m09',    'm10',    'm11',   'm12v',   'm12i',  'm12q', 'm13', 'm14']


    #sdirs=['m09_dw','m10_dw','m11_hr','m11_lr','hires_bh','m12_mr','m13_ics']
    #snums=[ 339,      440   ,  440  ,   440,    441    ,    441    ,   66  ]
    #bspecs=['-1',      '-1'  ,   '11',   '11',   '12',       '12'    ,   '13']
    #ymins=[ 1.e-4   ,  0.00005,   0.008,  0.03,   0.1   ,      0.3   ,   0.3]
    #ymaxs=[ 0.02,        0.02,   0.8,      3.,    9.   ,    40.     ,   40.]
    #rcuts=[ 15.,          15.,    15.,    5.,      20.  ,    5.     ,    50.]        

    #sdirs=['m11_hr']
    #snums=[440]
    #bspecs=['11']
    #ymins=[0.008]
    #ymaxs=[0.8]
    #rcuts=[15.]

    #sdirs=['m13_mr']
    #snums=[437]
    #bspecs=['13']
    #ymins=[0.3]
    #ymaxs=[100.]
    #rcuts=[15.]

    for i,slbl,sdir,snum,bspec,ymin,ymax,rcut in \
      zip(range(np.array(sdirs).size),slbls,sdirs,snums,bspecs,ymins,ymaxs,rcuts):
        if (i>-1):
            print(sdir)
            
            z=np.arange(0.,20.,0.1); sfr=10.**(np.random.rand(z.size));
            z,sfr = sfr_time_snap( sdir0+sdir, snum, rcut=rcut, dt=0.001, cen=[x000[i],y000[i],z000[i]])
            mcorr_mloss = 1.0 #0.7 + 0.3/(1.+z)
            if(slbl=='m12q'): mcorr_mloss = 0.7
            sfr /= mcorr_mloss
            xx = z
            xx = 1.+z
            #pylab.subplot(2,3,i+1)
            ii=i+1; 
            if(slbl=='m13'): ii+=1;
            #pylab.subplot(2,4,i+1)
            pylab.subplot(3,3,ii)

            lpos=(9.8,22.)
            if(slbl=='m09'): lpos=(1.1,1.e-2)
            if(slbl=='m10'): lpos=(9.8,1.e-2)
            if(slbl=='m11'): lpos=(9.8,1.7)
            if(slbl=='m12v'): lpos=(8.9,5.2)
            if(slbl=='m12q'): lpos=(8.9,22.)
            if(slbl=='m13'): lpos=(9.8,22.)
            if(slbl=='m13'): lpos=(9.8,49.2)
            lpos=(9.4,ymax*10.**(-0.09*np.log10(ymax/ymin)))
            pylab.text(lpos[0],lpos[1],slbl,color='k',fontsize=18)
            
            zmax = 8.
            zmax = 0.8
            zmax = 12.
            #pylab.axis([0.,zmax,ymin,ymax])
            pylab.axis([1.+0.,1.+zmax,ymin,ymax])
            #pylab.xscale('linear'); pylab.yscale('log')
            pylab.xscale('log'); pylab.yscale('log')
            ztick=1.+np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
            pylab.xticks(ztick,['0','1','2','3','4','5','6','7','8','','10','',''])
            #if(i==0)|(i==3):
            if(i==0)|(i==3)|(i==6):
                pylab.ylabel(r'${\rm SFR}\ \ \dot{M}_{\ast}\ \  [{\rm M_{\odot}\,yr^{-1}}]$')
            else:
                pylab.ylabel(r' ')
            #if(i>=3):
            if(i==3)|(i==5)|(i==6):
                pylab.xlabel(r'Redshift $z$')
            else:
                pylab.xlabel(r' ')
                
            #pylab.plot(xx,sfr,'b:')
            sfr_smoo_1=util.smooth(sfr, window_len=101, window='flat')
            #sfr_smoo_1=1.*sfr
            pylab.plot(xx,sfr_smoo_1,'k-')
            len0=101
            #len0=11
            len0=1001
            if (sfr.size<len0): len0=sfr.size-5
            sfr_smoo_2=util.smooth(sfr, window_len=len0, window='flat')
            pylab.plot(xx,sfr_smoo_2,'r--')
            #np.savetxt(sdir+'_sfr_maingal.txt',[,])
            asciitable.write({'z': xx, 'sfr': sfr_smoo_1}, sdir+'_sfr_maingal.txt', names=['z', 'sfr'])

            if (bspec != '-1'):
                bname=bdir0+'individual_tracks_corrected_'+bspec+'.0.dat'
                v=np.loadtxt(bname)
                a=v[:,0]; sfr_p=v[:,1]+v[:,2]; sfr_m=v[:,1]-v[:,3]; z=1./a-1.
                x=1.+z
                pylab.plot(x,sfr_p,'c-')
                pylab.plot(x,sfr_m,'c-')
                
            moster_dir='/Users/phopkins/Documents/work/plots/zooms/moster_sfhs/'
            bspec_m='-1'
            if(slbl=='m11'): bspec_m='11'
            if(slbl=='m12v'): bspec_m='12'
            if(slbl=='m12q'): bspec_m='12'
            if(slbl=='m12i'): bspec_m='12'
            if(slbl=='m13'): bspec_m='13'
            if(slbl=='m14'): bspec_m='14'
            if (bspec_m != '-1'):
                bname=moster_dir+'m'+bspec_m+'_plus'
                v=np.loadtxt(bname); oneplusz=np.array(v[:,0]); sfr=10.**(np.array(v[:,1]));
                pylab.plot(oneplusz,sfr,'m-')
                bname=moster_dir+'m'+bspec_m+'_minus'
                v=np.loadtxt(bname); oneplusz=np.array(v[:,0]); sfr=10.**(np.array(v[:,1]));
                pylab.plot(oneplusz,sfr,'m-')

            pylab.subplots_adjust(left=(18./25.4)/8., 
                bottom = (28./25.4)/12., 
                right = 1 - (2./25.4)/8., 
                top = 1 - (5./25.4)/12.,
                hspace = 0.121, wspace=0.112)

    pylab.savefig("zoom_sfh_z0_vsobs.pdf",transparent=True,bbox_inches='tight',pad_inches=0)



def oplot_beh():
    bdir0='/Users/phopkins/Documents/work/plots/python/behroozi_sfh/sfh_individual/'
    bspec='12'
    bname=bdir0+'individual_tracks_corrected_'+bspec+'.0.dat'
    v=np.loadtxt(bname)
    a=v[:,0]; sfr_p=v[:,1]+v[:,2]; sfr_m=v[:,1]-v[:,3]; z=1./a-1.
    x=z
    x=np.log(1.+z)
    pylab.plot(x,sfr_p,'b-')
    pylab.plot(x,sfr_m,'b-')


def sfh_z0_vs_subgrid():
    pylab.close('all'); 
    plot.figure(1,figsize=(8.,6.))
    charsi=22.
    matplotlib.rcParams.update({'font.size':18})

    sdir0='/Users/phopkins/Documents/work/plots/zooms/'

    zr=np.array([0.,8.]); pylab.xscale('linear');
    zr=np.array([0.,12.]); pylab.xscale('linear');
    zr=1.+zr; pylab.xscale('log');
    #pylab.axis(np.concatenate((zr,np.array([0.,20.])))); pylab.yscale('linear');
    #pylab.axis(np.concatenate((zr,np.array([0.1,50.])))); pylab.yscale('log');
    pylab.axis(np.concatenate((zr,np.array([0.02,50.])))); pylab.yscale('log');
   
    pylab.xlabel(r'Redshift $z$',fontsize=charsi)
    pylab.ylabel(r'${\rm SFR}\ \ \dot{M}_{\ast}\ \  [{\rm M_{\odot}\,yr^{-1}}]$',fontsize=charsi)
    pylab.xticks([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.],\
        ['0','1','2','3','4','5','6','7','8','9','10','',''])
    pylab.yticks([0.1,1.,10.],['0.1','1','10'])

    sdirs=['hires_bh','z10n192_B1hr','z10n192_B1hr_cwef1en25','z10n192_B1hr_cwef2en50']
    snums=[441,         441,            441,                    441                  ]

    sdirs=['B1_Aug27_egylim4','B1_Sep10_mhd','B1_Aug27_egylim','B1_Aug8']
    snums=[440,         440,            440,                    440                 ]

    sdirs=['B1_Aug27_egylim','z10n192_B1hr','z10n192_B1hr_cwef1en25','z10n192_B1hr_cwef2en50']
    snums=[440,         441,            441,                    441                  ]


    sdirs=['B1_11_Jan2014','z10n192_B1hr','z10n192_B1hr_cwef1en25','z10n192_B1hr_cwef2en50']
    snums=[440,         441,            441,                    441                  ]


    plotstyles=[ 'k-',   'r--',          'b-.',                 'g:'                ] 
    labels=[r'Explicit Feedback',r'No Feedback',r'Sub-Grid Wind 1',r'Sub-Grid Wind 2']
    for sdir,snum,plotstyle,labeler in zip(sdirs,snums,plotstyles,labels):
        print('crunching sdir=',sdir)
        z,sfr = sfr_time_snap( sdir0+sdir, snum, rcut=5., dt=0.02 )
        pylab.plot(1.+z,sfr,plotstyle,linewidth=2.,label=labeler)

    pylab.legend(loc='upper left',fontsize=10,frameon=False,labelspacing=0.1,numpoints=1)
    pylab.savefig("zoom_sfh_z0_subgrid.pdf",transparent=True,bbox_inches='tight',pad_inches=0)


def sfr_time_snap( sdir, snum, rcut=10., dt=0.01 , cen=[0.,0.,0.]):
    #cen = calculate_zoom_center(sdir,snum);
    P=gadget.readsnap(sdir,snum,4,cosmological=1,skip_bh=1);
    ok=(P['m']>0.)
    if(cen[0] != 0.):
        kc=1000.
        pos=P['p']; 
        x=pos[:,0]-cen[0]*kc; y=pos[:,1]-cen[1]*kc; z=pos[:,2]-cen[2]*kc;
        r=np.sqrt(x*x+y*y+z*z);
        ok=(r < rcut);
    a_form=P['age'][ok]; m=P['m'][ok];
    s=np.argsort(a_form); a_form=a_form[s]; m=np.cumsum(m[s]);
    dz=1.e-3; z_grid=np.arange(0.,30.,dz); a_grid=1./(1.+z_grid);
    #dz=1.e-2; z_grid=np.arange(0.,30.,dz); a_grid=1./(1.+z_grid);
    #dz=1.e-1; z_grid=np.arange(0.,30.,dz); a_grid=1./(1.+z_grid);
    t_grid=cosmo.lookback_time(z_grid).value;
    ok=(a_grid >= a_form[0]) & (a_grid <= a_form[-1]);    
    z_grid=z_grid[ok]; a_grid=a_grid[ok]; t_grid=t_grid[ok];
    t_grid_new = np.arange(np.max(t_grid),np.min(t_grid),-dt)
    a_grid_new = np.interp(t_grid_new, t_grid, a_grid)

    m_grid = np.interp(a_grid_new, a_form, m)
    N_grid = m_grid.size
    dm = np.concatenate((np.array([m_grid[0]]),m_grid[1:N_grid]-m_grid[0:N_grid-1]))
    dm_dt = dm / dt
    dm_dt *= (1.e10/1.e9); # units and mass loss
    z = 1./a_grid_new - 1.
    return z, dm_dt



def sfh_z0_daddi():
    pylab.close('all'); 
    sdir0='/Users/phopkins/Documents/work/plots/zooms/'
    bdir0='/Users/phopkins/Documents/work/plots/python/behroozi_sfh/sfh_individual/'

    plot.figure(1,figsize=(24.,6.))
    charsi=20.
    matplotlib.rcParams.update({'font.size':18})

    sdirs=['hires_bh','m12_mr','m13_ics','m14_tst']
    txtnm=['m11.5_sfr','m12.0_sfr','m13.0_sfr','m14.0_sfr']
    snums=[441    ,    441    ,   66, 275]
    bspecs=['12',       '12'    ,   '13', '14']
    ymins=[  0.1   ,      0.3   ,   0.3,  1.0]
    ymaxs=[  9.   ,    40.     ,   40., 1000.]
    rcuts=[ 20.  ,    5.     ,    50.,    100.]

    for i,sdir,snum,bspec,ymin,ymax,rcut in \
      zip(range(np.array(sdirs).size),sdirs,snums,bspecs,ymins,ymaxs,rcuts):
        if (i>-1):
            print(sdir)
            z,sfr = sfr_time_snap( sdir0+sdir, snum, rcut=15., dt=0.01 )
            mcorr_mloss = 0.7 + 0.3/(1.+z)
            sfr /= mcorr_mloss
            xx = z
            xx = 1.+z
            pylab.subplot(2,3,i+1)
            zmax = 18.
            pylab.axis([1.+0.,1.+zmax,ymin,ymax])
            pylab.xscale('log'); pylab.yscale('log')
            ztick=1.+np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
            pylab.xticks(ztick,['0','1','2','3','4','5','6','7','8','','10','',''])
            pylab.xlabel(r'Redshift $z$')
            pylab.ylabel(r'${\rm SFR}\ \ \dot{M}_{\ast}\ \  [{\rm M_{\odot}\,yr^{-1}}]$')
            pylab.plot(xx,sfr)
            q=np.zeros((z.size,2)); q[:,0]=z[:]; q[:,1]=sfr[:];
            np.savetxt(txtnm[i],q,fmt='%.4e')

    pylab.savefig("zoom_sfh_z0_daddi.pdf",transparent=True,bbox_inches='tight',pad_inches=0)
