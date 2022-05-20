from __future__ import division,print_function

## this file is just a pure convenience wrapper for the files in the visualization module

######################################################################

import visualization.contour_makepic

def contour_makepic( x, y, z, hsml, weight, 
        weight2=0, weight3=0, 
        xlen = 1, 
        pixels = 720, set_aspect_ratio = 1.0, 
        set_maxden = 1.0e-1, ## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
        set_dynrng = 1.0e4, 
        invert_colorscale = 0, 
        set_percent_maxden = 0, set_percent_minden = 0 ):
        
    return visualization.contour_makepic.contour_makepic(\
        x,y,z,hsml,weight,weight2=weight2,weight3=weight3,xlen=xlen,\
        pixels=pixels,set_aspect_ratio=set_aspect_ratio,set_maxden=set_maxden,set_dynrng=set_dynrng,\
        invert_colorscale=invert_colorscale,set_percent_maxden=set_percent_maxden,\
        set_percent_minden=set_percent_minden);
        
def simple_makepic( x, y, 
    xrange=[-1.,1.], yrange=0, weights=0, hsml=0, 
    dont_make_plots=0, color_temperature=0, temp_weights=0, 
    set_temp_max=0, set_temp_min=0, 
    pixels = 720, invert_colorscale = 0, 
    set_maxden = 1.0e-1, set_dynrng = 1.0e4, 
    set_percent_maxden = 0, set_percent_minden = 0 ):

    return visualization.contour_makepic.simple_makepic( x, y, 
        xrange=xrange, yrange=yrange, weights=weights, hsml=hsml, 
        dont_make_plots=dont_make_plots, color_temperature=color_temperature, temp_weights=temp_weights, 
        set_temp_max=set_temp_max, set_temp_min=set_temp_min, 
        pixels=pixels, invert_colorscale=invert_colorscale, 
        set_maxden=set_maxden, set_dynrng=set_dynrng, 
        set_percent_maxden=set_percent_maxden, set_percent_minden=set_percent_minden );
      
######################################################################

import visualization.get_attenuated_stellar_luminosities as vgsl

def get_attenuated_stellar_luminosities( BAND_IDS, star_pos, gas_pos, bh_pos, \
        stellar_age, stellar_metallicity, stellar_mass, \
        gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
        bh_luminosity, \
        xrange=0, yrange=0, zrange=0, \
        INCLUDE_BH=0, SKIP_ATTENUATION=0, 
        ADD_BASE_METALLICITY=0., ADD_BASE_AGE=0., 
        IMF_SALPETER=0, IMF_CHABRIER=1, \
        MIN_CELL_SIZE=0.01, OUTER_RANGE_OF_INT=1200., \
        SCATTERED_FRACTION=0.0, \
        REDDENING_SMC=0, REDDENING_LMC=0, REDDENING_MW=0, \
        AGN_MARCONI=0, AGN_HRH=1, AGN_RICHARDS=0, AGN_SDSS=0 ):

    return vgsl.get_attenuated_stellar_luminosities( \
        BAND_IDS, star_pos, gas_pos, bh_pos, \
        stellar_age, stellar_metallicity, stellar_mass, \
        gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
        bh_luminosity, \
        xrange=xrange, yrange=yrange, zrange=zrange, \
        INCLUDE_BH=INCLUDE_BH, SKIP_ATTENUATION=SKIP_ATTENUATION, 
        ADD_BASE_METALLICITY=ADD_BASE_METALLICITY, ADD_BASE_AGE=ADD_BASE_AGE, 
        IMF_SALPETER=IMF_SALPETER, IMF_CHABRIER=IMF_CHABRIER, \
        MIN_CELL_SIZE=MIN_CELL_SIZE, OUTER_RANGE_OF_INT=OUTER_RANGE_OF_INT, \
        SCATTERED_FRACTION=SCATTERED_FRACTION, \
        REDDENING_SMC=REDDENING_SMC, REDDENING_LMC=REDDENING_LMC, REDDENING_MW=REDDENING_MW, \
        AGN_MARCONI=AGN_MARCONI, AGN_HRH=AGN_HRH, AGN_RICHARDS=AGN_RICHARDS, AGN_SDSS=AGN_SDSS );

######################################################################

import visualization.make_threeband_image

def make_threeband_image( x, y, lums, hsml=0, xrange=0, yrange=0, \
    dont_make_image=0, maxden=0, dynrange=0, pixels=720, \
    color_scheme_nasa=1, color_scheme_sdss=0 ):

    return visualization.make_threeband_image.make_threeband_image( \
        x, y, lums, hsml=hsml, xrange=xrange, yrange=yrange, \
        dont_make_image=dont_make_image, maxden=maxden, dynrange=dynrange, pixels=pixels, \
        color_scheme_nasa=color_scheme_nasa, color_scheme_sdss=color_scheme_sdss );

######################################################################

import visualization.return_columns_to_sources as vrcts

def return_columns_to_sources( source_pos, gas_pos, \
    gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
    xrange=0, yrange=0, zrange=0, \
    MIN_CELL_SIZE=0.01, OUTER_RANGE_OF_INT=1200., \
    TRIM_PARTICLES=1 ):

    return vrcts.return_columns_to_sources( source_pos, gas_pos, \
        gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
        xrange=xrange, yrange=yrange, zrange=zrange, \
        MIN_CELL_SIZE=MIN_CELL_SIZE, OUTER_RANGE_OF_INT=OUTER_RANGE_OF_INT, \
        TRIM_PARTICLES=TRIM_PARTICLES );

######################################################################
      
import visualization.colors      

def rgb_to_hls(x):
    return visualization.colors.rgb_to_hls_v(x);
    
def hls_to_rgb(x):
    return visualization.colors.hls_to_rgb_v(x);

def hue_minus_convert(hue,lightness,saturation):
    return visualization.colors.hue_minus_convert(hue,lightness,saturation);

def invertreverse_rgb(R,G,B):
    return visualization.colors.invertreverse_rgb(R,G,B);

def temperature_map_color_index(mass_pic, temp, set_temp_max=0, set_temp_min=0, 
        huem100=0, invertreverse=0):
    return visualization.colors.temperature_map_color_index( mass_pic, temp, \
        set_temp_max=set_temp_max, set_temp_min=set_temp_min, huem100=huem100, invertreverse=invertreverse);

def load_my_custom_color_tables():
    return visualization.colors.load_my_custom_color_tables();

######################################################################
      