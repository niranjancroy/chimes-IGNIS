#niranjan/May 2022: tree code for a faster calculation of incedent stellar flux on the gas particles. The structure of the code is from PYTREEGRAC developed  by Mike Grudic (https://github.com/mikegrudic/pytreegrav.git) modified for flux calculations.

import numpy as np
from numpy import sqrt, empty, zeros, empty_like, zeros_like
from numba import njit, prange
from .octree import *
from .dynamic_tree import *


def FluxTarget_bruteforce(x_target, x_source, L_source):
    """Returns the flux at a set of target positions, due to a set of source particles.

    Arguments:
    x_target -- shape (N,3) array of positions where the flux is to be evaluate
    x_source -- shape (M,3) array of positions of star particles
    L_source -- shape (M,) array of source luminosity

    Returns:
    shape (N,3) array of stellar flux
    """
    flux = zeros_like(x_target[:,0])
    for i in prange(x_target.shape[0]):
        dx = zeros(3)
        for j in range(x_source.shape[0]):
            r2 = 0
            for k in range(3):
                dx[k] = x_source[j,k]-x_target[i,k]
                r2 += dx[k]*dx[k]
            if r2 == 0: continue # no force if at the origin
            r = sqrt(r2)
     
            flux[i] += L_source[j]*(1/r2)*(1/(4*np.pi))

    return flux

FluxTarget_bruteforce_parallel = njit(FluxTarget_bruteforce,fastmath=True,parallel=True)
FluxTarget_bruteforce = njit(FluxTarget_bruteforce,fastmath=True)

def ConstructTree(pos,L,softening=None,quadrupole=False,vel=None):
    """Builds a tree containing particle data, for subsequent potential/field evaluation

    Parameters
    ----------
    pos: array_like
        shape (N,3) array of particle positions
    m: array_like
        shape (N,) array of particle luminosities
    softening: array_like or None, optional
        shape (N,) array of particle softening lengths - these give the radius of compact support of the M4 cubic spline mass distribution of each particle
    quadrupole: bool, optional
        Whether to store quadrupole moments (default False)
    vel: bool, optional
        Whether to store node velocities in the tree (default False)

    Returns
    -------
    tree: octree
        Octree instance built from particle data
    """
    if softening is None: softening = zeros_like(L)
    if not (np.all(np.isfinite(pos)) and np.all(np.isfinite(L)) and np.all(np.isfinite(softening))):
        print("Invalid input detected - aborting treebuild to avoid going into an infinite loop!")
        raise
    if vel is None:
        return Octree(pos,L,softening,quadrupole=quadrupole)
    else:
        return DynamicOctree(pos, L,softening,vel, quadrupole=quadrupole)




@njit(fastmath=True)
def FluxWalk(pos, tree, softening=0, no=-1, theta=0.7):
    """Returns the stellar flux at position x by performing the Barnes-Hut treewalk using the provided octree instance
    Arguments:
    pos - (3,) array containing position of interest
    tree - octree instance storing the tree structure    
    Keyword arguments:
    softening - softening radius of the particle at which the force is being evaluated - we use the greater of the target and source softenings when evaluating the softened potential
    no - index of the top-level node whose field is being summed - defaults to the global top-level node, can use a subnode in principle for e.g. parallelization
    theta - cell opening angle used to control force accuracy; smaller is slower (runtime ~ theta^-3) but more accurate. (default 0.7 gives ~1% accuracy)
    """
    if no < 0: no = tree.NumParticles # we default to the top-level node index
    flux = 0 #zeros(1,dtype=np.float64)
    dx = np.empty(3,dtype=np.float64)

    h = max(tree.Softenings[no],softening) #this is defined just for the purpose of tree calculations, we'll have zero softening always

    while no > -1: # loop until we get to the end of the tree
        r2 = 0
        for k in range(3):
           dx[k] = tree.Coordinates[no,k] - pos[k]
           r2 += dx[k]*dx[k]

        r = sqrt(r2)
        sum_field = False

        if no < tree.NumParticles: # if we're looking at a leaf/particle
            if r > 0:  # no self-force
              # use point mass force
                fac = tree.Luminosities[no]/(r2) #niranjan
                sum_field = True
            no = tree.NextBranch[no]
        elif r > max(tree.Sizes[no]/theta + tree.Deltas[no], h+tree.Sizes[no]*0.6+tree.Deltas[no]): # if we satisfy the criteria for accepting the monopole            
            fac = tree.Luminosities[no]/(r2)
            sum_field = True
            no = tree.NextBranch[no] # go to the next branch in the tree
        else: # open the node
            no = tree.FirstSubnode[no]
            continue

        if sum_field: # OK, we have fac for this element and can now sum the flux
            flux += fac*(1/(4*np.pi)) 
  
    return flux


def FluxTarget_treecal(pos_target, softening_target, tree, theta=0.7):

      if softening_target is None: softening_target = zeros(pos_target.shape[0])

      #result = np.zeros(pos_target.shape[0], dtype=np.float64)
      result = zeros_like(pos_target[:,0])

      for i in prange(pos_target.shape[0]): 
            result[i] = FluxWalk(pos_target[i], tree, softening=softening_target[i], theta=theta)
   
      return result

FluxTarget_treecal_parallel = njit(FluxTarget_treecal,fastmath=True,parallel=True) 
FluxTarget_treecal = njit(FluxTarget_treecal,fastmath=True)


def FluxTarget_tree(pos_target, pos_source, L_source, theta=.7, softening_source=None, tree=None, return_tree=False, parallel=False):
    """Stellar flux calculation for general N+M body case

    Returns the flux for a set of M star particles with positions x_source and luminosities L_source, at the positions of a set of N particles that need not be the same.

    Parameters
    ----------
    pos_target: array_like
        shape (N,3) array of target particle positions where you want to know the flux
    pos_source: array_like
        shape (M,3) array of source particle positions (positions of particles sourcing the flux)
    L_source: array_like
        shape (M,) array of source particle luminosities
    softening_target: array_like or None, optional
        shape (N,) array of target particle softening radii - these give the radius of compact support of the M4 cubic spline flux distribution
    softening_source: array_like or None, optional
        shape (M,) array of source particle radii - these give the radius of compact support of the M4 cubic spline flux distribution
    theta: float, optional
        cell opening angle used to control force accuracy; smaller is slower (runtime ~ theta^-3) but more accurate. (default 0.7, gives ~1% accuracy)
    parallel: bool, optional
        If True, will parallelize the force summation over all available cores. (default False)
    tree: Octree, optional
        optional pre-generated Octree: this can contain any set of particles, not necessarily the target particles at pos (default None)
    return_tree: bool, optional
        return the tree used for future use (default False)

    Returns
    -------
    phi: array_like
        shape (N,1) array of flux at the target positions
    """

    quadrupole = False

    if softening_source is None and pos_source is not None: softening_source = np.zeros(len(pos_source))
    #if tree is None: tree = ConstructTree(np.float64(pos_source),np.float64(L_source), np.float64(softening_source), quadrupole=quadrupole) # build the tree if needed
    if tree is None: tree = ConstructTree(pos_source,L_source,softening_source, quadrupole=quadrupole) #niranjan: removing float64 as it's probably causing segfault due to high memory usage
    if parallel:
        flux = FluxTarget_treecal_parallel(pos_target, softening_target=None, tree=tree,theta=theta)
    else:
        flux = FluxTarget_treecal(pos_target, softening_target=None, tree=tree,theta=theta)

    return flux
