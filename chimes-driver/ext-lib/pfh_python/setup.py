import numpy as np
import matplotlib
#matplotlib.use('Agg') ## this calls matplotlib without a X-windows GUI
import matplotlib.pylab as pylab
import h5py as h5py
import os.path
import math
import scipy.interpolate as interpolate
import scipy.optimize as optimize
from astropy.cosmology import WMAP9 as cosmo
import gadget
import pfh_utils as util
import fit_functions as fitfun

fsize=20
pylab.rc("font", size=fsize, family='serif', weight='normal')
pylab.rc('text',usetex=True)
pylab.rc('font',**{'family':'serif','serif':['Computer Modern']})
pylab.rc("axes", labelsize=fsize, titlesize=1.5*fsize)
pylab.rc("xtick", labelsize=fsize)
pylab.rc("ytick", labelsize=fsize)
pylab.rc("legend", fontsize=fsize)
pylab.rc("axes", linewidth=1.0)
pylab.rc("lines", markeredgewidth=1)
pylab.tick_params('both', length=0.5*fsize, width=1, which='major')
pylab.tick_params('both', length=0.25*fsize, width=1, which='minor')
