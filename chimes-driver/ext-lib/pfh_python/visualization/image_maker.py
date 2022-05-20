import numpy as np
import matplotlib.pyplot as plot
from subprocess import call
import pfh_utils as util
from gizmopy.load_fire_snap import *
from gizmopy.quicklook import *
import visualization.colors as viscolors
import visualization.make_threeband_image as makethreepic
import visualization.contour_makepic as cmakepic
import visualization.raytrace_projection as rayproj



def edgeon_faceon_projection(sdir, snum, centering='', field_of_view=30., 
    edgeon=False, faceon=True,**kwargs):
    """
    routine to quickly generate edge-on or face-on image projections (single-image) 

    Args:
      sdir (str):  directory with snapshots (e.g. 'output' or 'm12i_res7100_md/output')
      snum (int):  snapshot number to visualize
      centering (string or len 3 array):  how to center the rotation & image.  see 
                                          image_maker docstring for a description.
      field_of_view (float):  field of view of the image
      edgeon (bool):  whether to make an edgeon image (otherwise, makes a faceon image)
      faceon (book):  whether to make a faceon image.   overridden by edgeon though (i.e.
                      edgeon == True and faceon == True yeilds an edgeon image)

      **kwargs: any valid kwargs for image_maker

    """

    # only load up the xyz particles for the primary ptype if we're doing COM centering
    if centering == 'com':
        image_maker_defaults = util.get_default_args(image_maker.image_maker)
        image_key = kwargs.get('image_key', image_maker_defaults['image_key'])
        if 'star' in image_key:
            ptype = image_maker_defaults['ptypes_stellar'][0]
        else:
            ptype = 0
        xyz_primary = load_fire_snap('Coordinates', ptype, sdir, snum, **kwargs)
    else:
        xyz_primary = None

    # parse the centering string / array into an array, possibly w/ estimate_zoom_center
    center = get_center(sdir, snum, xyz_primary, centering, **kwargs)
    x_vec,y_vec,z_vec = get_angmom_projection_vectors(sdir,snum,rcut=field_of_view/2.,center=center)
    proj=[x_vec,y_vec,z_vec]; imprefix='image_faceon'
    if(edgeon==True or faceon==False): 
        proj=[x_vec,z_vec,y_vec]; imprefix='image_edgeon'
    image_maker(sdir,snum,centering=center,field_of_view=field_of_view,
        projection_matrix=proj,image_name_prefix=imprefix,**kwargs)
    return


def image_maker(sdir, snum, image_key='star', 
    dynamic_range=1.e4, saturation_level=-0.0001,
    centering='', projection='',
    camera_opening_angle=90., field_of_view=30., 
    axis_ratio=1., euler_angles=[0,0,0], 
    projection_matrix=np.array([[1,0,0],[0,1,0],[0,0,1]]),
    pixels=1024, image_name_prefix='image', 
    output_directory='./', labels='scale_redshift',
    layer_cmap='heat_purple',layer_alpha=0.3, 
    colorscheme='NASA', disable_lighting=False,
    bands=np.zeros(0),kernel_width=np.zeros(0),
    ptypes_stellar=[4], depth_clipping_factor=10., 
    h_rescale_factor=3, h_threshold_factor=0.01,
    dust_to_gas_ratio_rescale=1., DesNgb=32, 
    preloaded_particle_data=None, **kwargs):
    '''
    Make ray-traced 3-band images of simulations.
    
    Core inputs::
      *sdir: parent directory containing the snapshots (string).  typically 'output'
      
      *snum: number (int) of snapshot (will automatically join multi-part snapshots)

      *image_key='string': tells the code what type of image to actually make
        'star' (or 'stellar') in string: then it will be a 3-band image of 
          starlight, attenuated by gas/dust, in the specified image bands
        'star_X' will create a multi-layer image, where the stellar image has an 
          additional image of gas property 'X' overlaid on it as a semi-transparent
          layer ('X' must be the name of a quantity which can be extracted with 
          'load_fire_snap', and will be represented as a single colormap with 
          increasing value: e.g. X = 'SpecificStarFormationRate', or 'Temperature' or 
          luminosities like 'CO' or 'Xray' or 'Halpha', or metallicity 'Z', etc.
        if 'star' does not appear, then the image will be a three-band image of the 
          gas alone, divided into three bands according to the specified gas property 
          and the values of the 'bands' parameter. default is Temperature, 
          i.e. a three-band image of the gas showing iso-density maps of gas around 
          three temperature ranges. but by setting image_key='Z' or 'Density' or 
          any other value readable by 'load_fire_snap', for example, you 
          can create a three-band map of metallicity around the desired iso-values
        examples: 'gas_Temperature' = 3-band gas temperature map
                  'star' = 'normal' HST-like stellar 3-band map
                  'star_Temperature' = stellar map with overlaid single-channel temp map

      *bands=[A,B,C] - define the bands used for the three-band images. 
        for starlight ('star' in image_key), these should be three bandpass names of the 
          allowed list, e.g. bands=['u','g','r'] (the default if this is not set).
          (allowed = 'Bol','U','B','V','R','I','J','H','K','u','g','r','i','z')
        for gas, these are the values of the band-centers for the scalar quantity the 
           gas map will represent. e.g. if image_key='gas_Temperature' (temperature plot),
           bands=[300., 2.0e4, 3.0e5] is the default, which are the temperatures 
           (in K) defining the 'cold', 'warm', and 'hot' phase centroids -- you can set 
           these to whatever values you like to best capture the phases of interest. 
           for other quantities, e.g. 'Density', you should set these to appropriate 
           values to separate the three bands/phases of interest
           
      *saturation_level - define the maximum surface brightness to be represented in the 
        image. behavior depends on whether this is >0 or <0 (positive or negative)
        - if set to a -positive- value: the code will assume this is a 
          specific, fixed, physical value, and use that as the maximum of the surface 
          brightness scale (so be sure to match this to your code units, e.g. for our 
          FIRE defaults this is in solar luminosities/kpc^2). if the whole image is 
          dark (saturated), this too high (low)
        - if set to a -negative- value: the code will assume this refers to the 
          fraction of pixels which should be saturated, and will dynamically scale 
          the intensity accordingly. default=-0.0001 means 0.01% of pixels are saturated
      
      *dynamic_range - dynamic range of the intensity scale between saturation, and 
        zero intensity (black). if too small, the image will be dark outside of the 
        brightest few pixels (if any). if too large, all the pixels will scrunch up
        against the saturation limit.  remember the intensity is log-scaled, though, 
        so this can be quite large (e.g. default=1e4, for showing faint 
        halo/low-surface-brightness features you may want to increase to 1e6, etc)


    Camera and centering options::
      *centering='string': determine the centering scheme for making the image. 
        default (not set) will be to center using the iterative center-finding scheme to 
          try and center the image on the center of the most massive galaxy/stellar object.
        set to a vector (e.g. centering=[1.,2.,3.]) to specify the center in 
          the same code units as the positions being read (forces a specific center)
        set to 'com' to center on the center-of-mass of the primary element type 
          (gas or stars) being imaged
        set to 'bh' to center on the median position of black hole particle[s] in the box
      
      *projection='string': determine the kind of projection scheme for the image.
        if 'camera' is not included in the string, it will be a flat projection (box is a 
          rectangular column, observer at infinite distance). 
        'camera': use a camera (with camera_opening_angle below), 
          pointing in the z-axis direction (this can be rotated as desired, below), 
          centered at the center point defined above
        'camera_tocenter': if 'tocenter' or 'to_center' included with camera, 
          a different convention is used where the camera is shifted along the 
          z-axis to that it points -at- the position defined by the centering above
      
      *camera_opening_angle=int: for camera projections, this sets the opening angle 
        of the lens (in degrees). default 90. wider gives more 'fisheye/all-sky' effects

      *field_of_view=float: set approximate field-of-view (in code units, e.g. kpc), 
        for images with fixed field ('flat' projection or 'camera_tocenter' pointed at center, 
        where this is the desired size at the central location (determines camera distance)

      *axis_ratio=float: specify the image-plane axis ratio =height/width (default=1)

      *euler_angles=[theta,psi,phi]=[pitch,roll,yaw]: Euler angles for pitch-yaw-roll rotations
        (default is [0,0,0], but use to rotate object about its center or rotate camera)
      
      *projection_matrix: arbitrary custom projection matrix given by
        [[x_axis_unit_vector],[y_axis_unit_vector],[z_axis_unit_vector]]
        (useful for custom camera systems, or easily getting exact face/edge-on images)

      
    Image style and file i/o options::
      *pixels=float: resolution of image in pixels on a side (will be pixels X pixels).
      
      *image_name_prefix='string': set default prefix for image filenames (default='image')
       (which use convention '[prefix]_s[snapshot number]_t[eulerangles pitch-yaw-roll]_Ngb[DesNgb]_[keyword for image].pdf'
    
      *output_directory='string': specify output folder for images. defaults to local
      
      *labels='string': specify if labels should be included by adding various terms to 
        this string. if it includes 'scale' a scale bar will be included (assumes kpc units)
        if it includes 'time' or 'Gyr' a time label (physical time units assumed, so dont 
        use this for cosmological simulations) will be included. if it includes 
        'z' or 'redshift' a redshift label will be included. if 'none' or '', no labels.
    
      *layer_cmap=colormap_name: give color map for extra semi-transparent image layers

      *layer_alpha=float: set alpha (0-1, 0=transparent, 1=opaque) for extra layers
      
      *colorscheme: set to 'SDSS' to use SDSS-style colors for stellar images, otherwise NASA 
      
      *disable_lighting: set True to force lighting effects always off (otherwise used in gas images)

      *kernel_width=[a,b,c] - (optionally) used with 'bands' if bands is set to numerical 
        values (not photometric bandpass names), in which case it defines the width 
        (in log_10-space) of the bandpass central (Gaussian) filter. default=[0.5,0.5,0.5]
    
      *ptypes_stellar - specify which particle types should be included as stars for 
        stellar images. default=[4], but often in non-cosmological simulations this 
        should be [2,3,4]


    For movie making/other specialized applications:
      *preloaded_particle_data:
        pass in a dictionary of particle data instead of loading it, e.g. if interpolating.
        items should be named as 'xyz' and 'xyz_g' for positions, but value that would be 
        passed to load_snapshot_vals (e.g. 'SmoothingLength') otherwise.  will 
        NB -- particle data should not be processed in any way -- center shift etc will 
        still be applied

      **kwargs:  passed to generate_directory_and_filename, load_fire_snap, and 
                 load_and_center_xyz_coordinates

    '''

    # setup the figure environment, create the working directory, define filenames
    image_key,effects_key,gasvalue_to_image,bands,kernel_width = parse_gasvals_bands_toplot(image_key,
        bands,kernel_width,disable_lighting=disable_lighting)

    filename = generate_directory_and_filename(sdir,snum,euler_angles=euler_angles,field_of_view=field_of_view,
        output_directory=output_directory,image_key=image_key,image_name_prefix=image_name_prefix, DesNgb=DesNgb,
        **kwargs)

    fig_axis = generate_figure_and_axes(axis_ratio=axis_ratio, pixels=pixels)

    # load xyz coordinates, center (using specified method), rotate and project
    if preloaded_particle_data is not None:
        time = preloaded_particle_data.pop('Time')
    else:
        time = load_fire_snap('Time',0,sdir,snum,**kwargs)
    m_g=0
    h_g=0
    z_g=0
    t_g=0

    ptypes=[0]
    if('star' in image_key): 
        ptypes=ptypes_stellar; # ptypes for main centering
    elif 'dark' in image_key:
        ptypes=[1]
    rot_matrix = np.dot(return_rotation_matrix(euler_angles),projection_matrix)

    xyz, xyz_g = load_and_center_xyz_coordinates(ptypes,sdir,snum,rot_matrix,
        centering=centering,image_key=image_key,preloaded_particle_data=
        preloaded_particle_data,**kwargs)  

    # compute this in case we print an error message later....
    mean_xyz = np.mean(xyz, axis=0)

    # project into camera and build clipping masks to keep only needed data
    xyz, xyz_g, ok, ok_g, xr, yr, zr = build_mask_and_camera(xyz,xyz_g,
        projection=projection,field_of_view=field_of_view,
        depth_clipping_factor=depth_clipping_factor,
        camera_opening_angle=camera_opening_angle,axis_ratio=axis_ratio)        

    if xyz.size == 0:
        print("!! no particles within field of view.  check your centering.")
        print("before the recenter,\n\tmean(xyz) = {:.2f}, {:.2f}, {:.2f}".format(*mean_xyz))
        print("with a centering method of")
        if type(centering) == str:
            print("\n\tcentering = {}".format(centering))
        else:
            print("\n\tcenter = {:.2f}, {:.2f}, {:.2f}".format(*centering))
        return None, None

    # load (clipped) data blocks needed for actual image generation    
    m,h,z,t = load_snapshot_brick(ptypes,sdir,snum,particle_mask=ok,
      gasvalue_to_image=gasvalue_to_image, DesNgb=DesNgb, 
      preloaded_particle_data=preloaded_particle_data)

    if('star' in image_key): 
        # also load the data for the gas to do attenuation
        m_g,h_g,z_g,t_g=load_snapshot_brick([0],sdir,snum,particle_mask=ok_g,
          gasvalue_to_image=gasvalue_to_image,preloaded_particle_data=preloaded_particle_data)
        
        # arbitrary re-scaling allowed by advanced options (occasionally useful for image tweaking)    
        z_g*=dust_to_gas_ratio_rescale
        h*=h_rescale_factor
        
        hmax=h_threshold_factor*np.abs(xr[1]-xr[0]); steepness=1.; h=1./(1./h**steepness + 1./hmax**steepness); # steepness = 2 gives much sharper cutoff, produces more 'point-like' images: looks more like point sources but less color fidelity        

    if('camera' in projection): m /= xyz[2]**2 + h**2 + zr[0]**2; # these need to be 'fluxes' (dimmer further)
    
    print("sending {:,} particles of the primary species to the ray tracer".format(m.size), flush=True)

    # ray-trace for image: these routines process the band-weights and call the ray-trace routines
    out_gas,out_u,out_g,out_r,massmap_xlayer,image_xlayer = process_raytrace_and_layercreation(image_key,
        effects_key,bands,kernel_width,xyz,m,t,z,h,xyz_g,m_g,t_g,z_g,h_g,xrange=xr,yrange=yr,zrange=zr,
        pixels=pixels, **kwargs)
    
    ## process the maps into actual images, add labels, scale-bars, clean up axes, etc
    image24, massmap = process_bandmaps_and_effects(fig_axis,image_key,effects_key,time,
        out_r,out_g,out_u,image_xlayer,colorscheme=colorscheme,labels=labels,
        dynamic_range=dynamic_range,xrange=xr,yrange=yr,layer_cmap=layer_cmap,
        layer_alpha=layer_alpha,pixels=pixels,saturation_level=saturation_level)

    # save the data
    plot.savefig(filename+'.pdf',dpi=pixels,bbox_inches='tight',pad_inches=0); 
    outfi=h5py.File(filename+'.hdf5','w')
    for name,dataset in zip(
        ['image24','massmap','image24_extralayer','massmap_extralayer','time','xrange','yrange',
         'image_key',  'effects','labels','colorscheme','extralayer_colormap','extralayer_alpha',
         'pixel_count','saturation_level','image_bands','image_projection','snapshot_directory',
         'snapshot_number','centering','camera_opening_angle','euler_angles','ProjectionMatrix','kernel_width'],
        [ image24 , massmap , image_xlayer       , massmap_xlayer     , time , xr     , yr     ,
          image_key ,effects_key, labels , colorscheme , layer_cmap          , layer_alpha      ,
          pixels      , saturation_level ,   bands     ,      projection  ,     sdir           ,
          snum            , centering , camera_opening_angle , euler_angles ,projection_matrix , kernel_width ]):
        outfi.create_dataset(name,data=dataset)
    outfi.close(); plot.close('all');
    print("-- saved {} --".format(filename))
    return


def process_raytrace_and_layercreation(image_key,effects_key,bands,kernel_width,
    xyz,m,t,z,h,xyz_g,m_g,t_g,z_g,h_g,xrange=[-1.,1.],yrange=[-1.,1.],zrange=[-1.,1.],
    pixels=750., openmp=True, **kwargs):
    ''' 
    call several subroutines to ray-trace the images. this is convenience entirely
    designed to localize import calls of other files

    n.b. -- kwargs don't do anything, but are included to allow any arguments

    '''
    image_xlayer=0; massmap_xlayer=0;
    if('star' in image_key):
        out_gas,out_u,out_g,out_r = rayproj.stellar_raytrace(bands,xyz[0],xyz[1],xyz[2],m,t,z,h,
            xyz_g[0],xyz_g[1],xyz_g[2],m_g,z_g,h_g,xrange=xrange,yrange=yrange,zrange=zrange,pixels=pixels,
            ADD_BASE_METALLICITY=2.e-8,ADD_BASE_AGE=3.e-4,IMF_SALPETER=False,IMF_CHABRIER=True, openmp=openmp)
        if('layer' in effects_key):
            massmap_xlayer,image_xlayer = cmakepic.simple_makepic(xyz_g[0],xyz_g[1],hsml=h_g,
                weights=m_g*t_g/np.sum(m_g*t_g),set_percent_maxden=0.9999,set_percent_minden=0.01,
                xrange=xrange,yrange=yrange,color_temperature=False,pixels=pixels,invert_colorscale=True)
    else: 
        out_gas,out_u,out_g,out_r = rayproj.gas_raytrace_temperature(bands,xyz[0],xyz[1],xyz[2],
            t,m,h,xrange=xrange,yrange=yrange,zrange=zrange,pixels=pixels,isosurfaces=True,
            kernel_width=kernel_width,add_temperature_weights=False,KAPPA_UNITS=2.0885*np.array([1.1,4.,2.]),
            openmp=openmp)
    return out_gas,out_u,out_g,out_r,massmap_xlayer,image_xlayer


def process_bandmaps_and_effects(fig_axis,image_key,effects_key,time,out_r,out_g,out_u,bonuslayer,
        colorscheme='NASA',labels='scale_time',xrange=[-1.,1.],yrange=[-1.,1.],
        dynamic_range=1.e4,saturation_level=-0.0001,layer_cmap='heat_purple',layer_alpha=0.3,**kwargs):
    ''' call several subroutines to process maps into an image, this is just convenience '''
    nasa_colors=True; sdss_colors=False; 
    if(colorscheme=='SDSS'): 
        nasa_colors=False; sdss_colors=True;
    maxden = saturation_level; ## default to assuming the maximum/saturation density is a fixed value
    if(saturation_level < 0.):
        ## but if we are setting the saturation level to a fractional (pixel-count) level, determine levels here

        # sum up the three band maps to get total intensities, and sort the pixels by intensity
        im_flat = np.sort(np.ndarray.flatten(out_r+out_g+out_u))

        if np.nanmax(im_flat) <= 0:
            print("max(r, g, u) = {:.2f}, {:.2f}, {:.2f}".format(np.nanmax(out_r), np.nanmax(out_g), np.nanmax(out_u)))
            print("max of flat image is zero; setting maxden to 1, but should be blank image I think")
            maxden = 1.0

        else:
          # chop down to only the portion of the image that's finite and > 0
          im_flat = im_flat[(im_flat>0.)&(np.isfinite(im_flat))]

          pixel_tot=out_r.size
          one_pix=1.00001/pixel_tot
          
          print('saturation_level == ',saturation_level,-one_pix,-(1.-one_pix))

          if(saturation_level > -one_pix): 
            saturation_level=-one_pix;
          if(saturation_level < -(1.-one_pix)): 
            saturation_level=--(1.-one_pix);

          # find the `saturation_level`th pixel in brightness and use it's brightness as maxden
          n0=np.round(np.abs(saturation_level)*im_flat.size).astype('int')
          if(n0<=0): 
            n0=0
          if(n0>=im_flat.size-1): 
            n0=im_flat.size-1
          maxden = im_flat[::-1][n0]
          if((maxden>=np.max(im_flat)) or (np.isnan(maxden))): 
            maxden=np.min(im_flat)
          if((maxden<=0) or (np.isnan(maxden))): 
            maxden=np.min(im_flat)     

    image24, massmap, im_max, im_min = makethreepic.make_threeband_image_process_bandmaps(out_r,out_g,out_u,
        maxden=maxden,dynrange=dynamic_range,color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors,**kwargs);
    
    print('saturation_threshold_in_code_units=',maxden,' dynamic_range=',dynamic_range)
    if('star' not in image_key): 
      image24 = makethreepic.layer_band_images(image24,massmap); 
    if('lighting' in effects_key): 
      image24 = overlay_lighting_layer(massmap,image24,maxden/dynamic_range,maxden)
    
    plot.imshow(image24,origin='lower',interpolation='bicubic'); # generates an image!

    if('layer' in effects_key): # overlay multi-layer images if desired
        viscolors.load_my_custom_color_tables();
        plot.imshow(bonuslayer,origin='lower',interpolation='bicubic',cmap=layer_cmap,alpha=layer_alpha)

    ## slap on the figure labels/scale bars
    if('scale' in labels or 'Scale' in labels): 
        overlay_scale_label(xrange,yrange,fig_axis)
    if('Gyr' in labels or 'Time' in labels or 'time' in labels): 
        overlay_time_label(time,fig_axis,label_redshift=False)    
    if('z' in labels or 'Redshift' in labels or 'redshift' in labels): 
        overlay_time_label(time,fig_axis,label_redshift=True)

    plot.subplots_adjust(left=0.0,bottom=0.0,right=1-0.0,top=1-0.0,hspace=0.0,wspace=0.0);
    fig_axis.axes.xaxis.set_ticklabels([]) #no tick names
    fig_axis.axes.yaxis.set_ticklabels([]) #no tick names
    fig_axis.axes.get_xaxis().set_ticks([]) #no ticks
    fig_axis.axes.get_yaxis().set_ticks([]) #no ticks
    plot.tight_layout(pad=0,w_pad=0,h_pad=0)

    return image24, massmap
    

def overlay_lighting_layer(massmap,image24,min,max,vmin_ratio=-1./6.,vmax_ratio=1./6.):
    ''' overlay lighting effects to give added depth and bring out features in some maps (takes playing around)'''
    light = viscolors.CustomLightSource(azdeg=0,altdeg=65)
    if (len(massmap.shape)>2): # three-band image, need to do more sophisticated treatment
        ## do some clipping to regulate the lighting:
        elevation=(massmap.sum(axis=2)-min)/(max-min); elevation[elevation<0.]=0.; elevation[elevation>1.]=1.;
        image24_lit = light.shade_rgb(image24,elevation*max,vmin=max*vmin_ratio,vmax=max*vmax_ratio)
    else: # single-band image, trivial shader here
        image24_lit = light.shade(image24, massmap)
    return image24_lit

      
def parse_gasvals_bands_toplot(image_key, bands, kernel_width, disable_lighting=False, **kwargs):
    ''' subroutine to set sensible defaults, and parse which types of quantities to plot '''
    if('Star' in image_key or 'stellar' in image_key or 'Stellar' in image_key): image_key+='star'
    gasval='Temperature'; bands=np.array(bands); kernel_width=np.array(kernel_width); 
    if('star' in image_key):
        if(bands.size <= 1): bands=np.array(['u','g','r'])
        band_dict={"Bol":0,"U":1,"B":2,"V":3,"R":4,"I":5,"J":6,"H":7,"K":8,"u":9,"g":10,"r":11,"i":12,"z":13}; bands_new=np.array([0,0,0]);
        for k in [0,1,2]: bands_new[k]=band_dict[bands[k]];
        bands=bands_new
    if(('star' not in image_key) or ('_' in image_key)):
        s0=image_key.split("_"); gasval=s0[len(s0)-1];
        if(gasval=='temp'): gasval='Temperature'
        if(gasval=='rho'): gasval='Density'
    if(('star' not in image_key) and (bands.size <= 1)):
        bands=np.array([300., 2.0e4, 3.0e5 ]) # default, like temperature
        if(gasval=='Density'): bands=np.array([10., 1. , 1.e-2])/406.25 # convert to nH in CGS, dense/WIM/CGM
    if(kernel_width.size <= 1): kernel_width=np.array([0.6,0.4,0.6])     
    effects_key=''
    if('star' in image_key and '_' in image_key): effects_key+='layer'
    if(('star' not in image_key) and (disable_lighting==False)): effects_key+='lighting'
    return image_key, effects_key, gasval, bands, kernel_width


def load_snapshot_vals(value,ptypes,sdir,snum,**kwargs):
    ''' subroutine to call load_fire_snap to obtain relevant values, 
    and automatically concatenate over a list of particle types if necessary '''
    matrix_key = False; ptypes = np.array(ptypes);
    for i, ptype in enumerate(ptypes):
        V=load_fire_snap(value,ptype,sdir,snum,**kwargs)
        if(i==0): 
            V0=1.*V; 
            if(np.array(V.shape).size > 1): matrix_key = True
        else:
            if(matrix_key): 
                V0=np.concatenate((V0,np.array([V])))
            else:
                V0=np.concatenate((V0,np.array(V)))
    return V0


def get_ptype_val_from_preloaded_pdata(preloaded_particle_data, ptypes, value, **kwargs):
  """
  subroutine to pop out data from a dictionary that's passed in
  """
  if len(ptypes) > 1:
    raise NotImplementedError("Only know how to handle a single particle type here, sorry")

  prefix = ''
  if ptypes[0] == 0:
    prefix = 'gas_'

  data = preloaded_particle_data.pop(prefix+value)

  if 'particle_mask' in kwargs:
    data = data.take(kwargs['particle_mask'], axis=0)
  return data


def load_snapshot_brick(ptypes,sdir,snum,gasvalue_to_image='Temperature',
    preloaded_particle_data=None, **kwargs):
    ''' subroutine to load and return the block of needed snapshot data '''
    if preloaded_particle_data is not None:
        h_sml = get_ptype_val_from_preloaded_pdata(preloaded_particle_data, ptypes, 'SmoothingLength', **kwargs)
        m = get_ptype_val_from_preloaded_pdata(preloaded_particle_data, ptypes, 'Masses', **kwargs)
        z_metallicity = get_ptype_val_from_preloaded_pdata(preloaded_particle_data, ptypes, 'Z', **kwargs)
        if(0 in ptypes):
            x = get_ptype_val_from_preloaded_pdata(preloaded_particle_data, ptypes, gasvalue_to_image, **kwargs)
            if('CO' in gasvalue_to_image or 'Xray' in gasvalue_to_image or 'Halpha' in gasvalue_to_image): x/=m; # make specific values, leads to correct weighting
        else:
            x = get_ptype_val_from_preloaded_pdata(preloaded_particle_data, ptypes, 'StellarAgeGyr', **kwargs)
    else:
        h_sml = load_snapshot_vals('SmoothingLength',ptypes,sdir,snum,**kwargs)
        m = load_snapshot_vals('Masses',ptypes,sdir,snum,**kwargs)
        z_metallicity = load_snapshot_vals('Z',ptypes,sdir,snum,**kwargs)
        if(0 in ptypes):
            x = load_snapshot_vals(gasvalue_to_image,ptypes,sdir,snum,**kwargs)
            if('CO' in gasvalue_to_image or 'Xray' in gasvalue_to_image or 'Halpha' in gasvalue_to_image): x/=m; # make specific values, leads to correct weighting
        else:
            x = load_snapshot_vals('StellarAgeGyr',ptypes,sdir,snum,**kwargs)
    return m, h_sml, z_metallicity, x


def get_center(sdir,snum,xyz_primary,centering='',**kwargs):
    """
    subroutine to determine center: keyword 'centering' defaults to 
    finding center, but if set to array it will use that, or if 'bh' 
    will use BH particle, or if 'com' it will use median position
    """
    cen=np.array(centering)
    if(cen.size>1): return cen; ## center is pre-determined for us
    if('bh' in cen): return np.median(load_fire_snap('Coordinates',5,sdir,snum,**kwargs),axis=0); ## center on median BH coordinates
    if('com' in cen): return np.median(xyz_primary,axis=0); ## center on median position

    print("computing center...", end='')
    cen, npart = estimate_zoom_center(sdir,snum,**kwargs) ## use fancy routine for centering
    if npart > 0:
        print(cen, flush=True)
        return cen
    else:
        print("!! warning:  couldn't center!  returning zeros")
        return np.zeros(3)


def load_and_center_xyz_coordinates(ptypes,sdir,snum,rot_matrix,centering='',image_key='star',
    preloaded_particle_data=None, **kwargs):
    ''' subroutine to extract, center, rotate the xyz coordinates as needed '''
    if preloaded_particle_data is not None:
      ## assume particle data has been passed in directly as a dictionary; pop out items as needed
      xyz = preloaded_particle_data.pop('xyz')
    else:
       # load coordinates of type of interest - need initial block
      xyz = load_snapshot_vals('Coordinates',ptypes,sdir,snum,**kwargs)

    center = get_center(sdir,snum,xyz,centering=centering,**kwargs) # sort out centering 
    xyz -= center; 
    xyz = xyz.transpose()
    xyz = np.dot(rot_matrix, xyz); # rotate matrix appropriately

    xyz_g = np.zeros(3)
    # handle gas if we're going stars as primary type
    if('star' in image_key): 
        if preloaded_particle_data is not None:
            xyz_g = preloaded_particle_data.pop('xyz_g')-center
        else:
            xyz_g=load_fire_snap('Coordinates',0,sdir,snum,**kwargs)-center
        xyz_g=xyz_g.transpose()
        xyz_g=np.dot(rot_matrix,xyz_g)
    return xyz, xyz_g


def build_mask_and_camera(xyz,xyz_g,projection='camera',field_of_view=1.,
        depth_clipping_factor=10.,camera_opening_angle=90.,axis_ratio=1.):
    '''without modification, the projection is 'flat': need to check for camera projection options.
    this routine projects to cameras of different styles, and builds masks for which elements to keep'''
    xr=0.5*field_of_view*np.array([-1.,1.]); 
    zr=depth_clipping_factor*xr; # assign 'comfort factor' to z-range to clip, to avoid too-distant field (warning: this clip will lead to pop-in)
    if('camera' in projection):
        xr=np.array([-1.,1.])*np.tan(camera_opening_angle*np.pi/360.); # define the window size in a useful angular coordinate
        z=xyz[2]; z_g=xyz_g[2]; c_dist=0.5*field_of_view/xr[1]; 
        if('tocenter' in projection or 'to_center' in projection): # normally the camera -is- the center, but if this is specified, camera -points- at center
            z+=c_dist; z_g+=c_dist; # move the particles 'along' the z-axis, effectively shifting the camera with respect to the center
        xyz[0]/=z; xyz_g[0]/=z_g; xyz[1]/=z; xyz_g[1]/=z_g; # correct to the angular coordinate
        zr=[0.1*c_dist,np.sqrt((2.*c_dist)**2 + (depth_clipping_factor*field_of_view)**2)] # comfortable factor to avoid dealing with much-too-distant objects
    yr=xr*axis_ratio
    safetyfac=1.15; xr*=safetyfac; yr*=safetyfac;

    ok=np.where((xyz[0] > xr[0])&(xyz[0] < xr[1])&(xyz[1] > yr[0])&(xyz[1] < yr[1])&(xyz[2] > zr[0])&(xyz[2] < zr[1]))[0]
    ok_g=np.where((xyz_g[0] > xr[0])&(xyz_g[0] < xr[1])&(xyz_g[1] > yr[0])&(xyz_g[1] < yr[1])&(xyz_g[2] > zr[0])&(xyz_g[2] < zr[1]))[0]
    xyz=xyz.take(ok,axis=1); 
    if(xyz_g[0].size > 1): xyz_g=xyz_g.take(ok_g,axis=1); # clip coordinates
    return xyz, xyz_g, ok, ok_g, xr/safetyfac, yr/safetyfac, zr


def return_perp_vectors(a):
    ''' procedure which will return, for a given input vector A_in, 
    the perpendicular unit vectors B_out and C_out which form perpendicular axes to A '''
    eps = 1.0e-10
    a = np.array(a,dtype='f');
    a /= np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    for i in range(len(a)):
        if (a[i]==0.): a[i]=eps;
        if (a[i]>=1.): a[i]=1.-eps;
        if (a[i]<=-1.): a[i]=-1.+eps;
    a /= np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    ax=a[0]; ay=a[1]; az=a[2];
    ## use a fixed rotation of the a-vector by 90 degrees:
    ## (this anchors the solution so it changes *continously* as a changes)
    t0=np.double(np.pi/2.e0);
    bx=0.*ax; by=np.cos(t0)*ay-np.sin(t0)*az; bz=np.sin(t0)*ay+np.cos(t0)*az;
    ## c-sign is degenerate even for 'well-chosen' a and b: gaurantee right-hand 
    ##  rule is obeyed by defining c as the cross product: a x b = c
    cx=(ay*bz-az*by); cy=-(ax*bz-az*bx); cz=(ax*by-ay*bx); 
    B_out=np.zeros(3); C_out=np.zeros(3);
    B_out[:]=[bx,by,bz]; C_out[:]=[cx,cy,cz];
    return B_out, C_out


def get_angmom_projection_vectors(sdir,snum,ptype=0,rcut=15.,center=[0.,0.,0.]):
    ''' quickly obtain the angular momentum vector and perpendicular vectors '''
    xyz,mask = center_and_clip(load_fire_snap('Coordinates',ptype,sdir,snum),center,rcut,1e9) # load positions
    m = load_fire_snap('Masses',ptype,sdir,snum).take(mask,axis=0) # load masses [or weights, more generically]
    vxyz = load_fire_snap('Velocities',ptype,sdir,snum).take(mask,axis=0) # load velocities
    vxyz = vxyz - np.median(vxyz,axis=0) # reset to median local velocity about r=0
    jvec=np.cross(vxyz,xyz); m_jvec=(m*(jvec.transpose())).transpose(); # compute j vector
    j_tot=np.sum(m_jvec,axis=0); z_vec=j_tot/np.sqrt(np.sum(j_tot*j_tot)); # this is the z-vector
    x_vec,y_vec = return_perp_vectors(z_vec)
    return x_vec, y_vec, z_vec    
        

def return_rotation_matrix(euler_angles):
   ''' give rotation matrix associated with euler angles 
        euler_angles=[theta,psi,phi]=[pitch-roll-yaw] in -degrees- '''
   tx,ty,tz = np.array(euler_angles) * np.pi/180.;
   Rx = np.array([[1,0,0], [0, np.cos(tx), -np.sin(tx)], [0, np.sin(tx), np.cos(tx)]])
   Ry = np.array([[np.cos(ty), 0, -np.sin(ty)], [0, 1, 0], [np.sin(ty), 0, np.cos(ty)]])
   Rz = np.array([[np.cos(tz), -np.sin(tz), 0], [np.sin(tz), np.cos(tz), 0], [0,0,1]])
   return np.dot(Rx, np.dot(Ry, Rz))


def generate_directory_and_filename(sdir, snum, euler_angles=[0,0,0], field_of_view=1., 
        DesNgb=32, output_directory='./', image_key='star', image_name_prefix='image', 
        format='.pdf', **kwargs):
    ''' generate names for outputs file[s] '''

    ## make sure output directory exists, and make any 
    ## required directories along the way (assuming permissions)
    os.makedirs(output_directory, exist_ok=True) 
    
    ## encode the snapshot number, field of view, and euler angles, 
    ## and number of neighbors in the output filename
    ss=snap_ext(snum,four_char=True); 
    fov=snap_ext(np.around(field_of_view).astype(int))

    ## but only include the euler angles if at least one is nonzero
    if (np.array(euler_angles) == np.zeros(3)).all():
      tt = ''
    else:
      tt='_t{}-{}-{}'.format(
        snap_ext(np.around(euler_angles[0]).astype(int)),
        snap_ext(np.around(euler_angles[1]).astype(int)),
        snap_ext(np.around(euler_angles[2]).astype(int)))

    fname_base = '{}/{}_s{}{}_fov{}_Ngb{}_{}'.format(
      output_directory, image_name_prefix, ss, tt, fov, DesNgb, image_key)

    return fname_base # return base name for figure

def generate_figure_and_axes(axis_ratio=1, pixels=1024):
    '''
    generage figure for output file
    '''
    plot.close('all'); # kill anything open
    fig=plot.figure(frameon=False,dpi=pixels,figsize=(1.,1.*axis_ratio)) # open new figure
    fig.set_size_inches(1.,1.*axis_ratio); # 1-inch (doesn't matter with PDF, all in dpi)
    ax_fig=plot.Axes(fig,[0.,0.,1.,1.*axis_ratio],clip_on=True) # max sure axes are set
    frame1=plot.gca() # identify axis constrols
    plot.xlabel(''); plot.ylabel('') # no axis labels
    for ax in [ax_fig,frame1]:
        ax.set_axis_off() # generally don't want to show the axes
        ax.axes.xaxis.set_ticklabels([]) # no tick names
        ax.axes.yaxis.set_ticklabels([]) # no tick names
        ax.axes.get_xaxis().set_ticks([]) # no ticks
        ax.axes.get_yaxis().set_ticks([]) # no ticks
    ax_fig=frame1 # assignment
    return ax_fig  # return axes object to use



def overlay_scale_label(xr,yr,figure_axis,color='w'):
    ''' function to overlay label for spatial scales, assume ~kpc units '''
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
    
    xlen=(xr[1]-xr[0]); ylen=(yr[1]-yr[0]); xoff = (0.25+0.02)*ddx / xlen; yoff = 0.025*(yr[1]-yr[0]) / ylen
    xr_new = np.array([xoff-dx/xlen*0.5,xoff+dx/xlen*0.5]); yr_new = np.array([yoff,yoff])
    plot.text(xoff,1.75*yoff,dxl,color=color,horizontalalignment='center',
        verticalalignment='baseline',transform=figure_axis.transAxes,fontsize=3)
    figure_axis.autoscale(enable=False,axis='both',tight=True)
    figure_axis.plot(xr_new,yr_new,color=color,linewidth=0.7,transform=figure_axis.transAxes)
    return
    
    
def overlay_time_label(time, figure_axis, label_redshift=True, color='w', n_sig_add=0, time_unit=' Gyr'):
    ''' function to overlay label for time, assumes ~Gyr scales '''
    time_to_use=time; prefix=''; suffix=time_unit; n_sig=2;
    if(label_redshift): time_to_use=1./time-1.; prefix='z='; suffix='' ## redshift
    if(time_to_use < 1.): n_sig += 1
    t_str = round_to_n( time_to_use, n_sig+n_sig_add )
    label_str = prefix+t_str+suffix
    plot.text(0.03,1.-0.025,label_str,color=color,horizontalalignment='left',
        verticalalignment='top',transform=figure_axis.transAxes,fontsize=3)
    return
    

def round_to_n(x, n):
    ''' Utility function used to round labels to significant figures for display purposes  '''
    if(n < 1): raise ValueError("number of significant digits must be >= 1")
    # show everything as floats (preference; can switch using code below to showing eN instead
    format = "%." +str(n-1) +"f"; as_string=format % x;
    return as_string
    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    if as_string[-3:] in ['+00', '+01', '+02', '+03','-01', '-02', '-03']:
        #then number is 'small', show this as a float
        format = "%." +str(n-1) +"f"; as_string=format % x
    return as_string

