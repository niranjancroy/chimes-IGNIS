.. Adding CHIMES
   Alexander Richings, 18th March 2020

.. _AddingChimes_label:

Adding CHIMES to a Hydro Code
=============================

In addition to running CHIMES as a stand-alone module with the CHIMES Driver python script, it can also be run on-the-fly within a hydrodynamics simulation to solve the non-equilibrium evolution of the chemical abundances alongside the hydrodynamics, gravity etc., and to use the resulting abundances to calculate the radiative cooling and heating rates and hence evolve the gas temperature. This involves replacing the cooling routines in the hydrodynamics code with calls to the CHIMES solver to integrate the chemical abundances and temperature for each gas particle every hydrodynamical time-step.

To add CHIMES to your own hydrodynamics code, you will need to add the source and header files from the ``src/`` directory in the ``chimes`` repository to your code. These will need to be added to the Makefile for your hydrodynamics code so that they will all be compiled and built together, and you will need the Sundials and HDF5 libraries. The ``Makefile_template`` in the ``chimes`` repository is used to build CHIMES as a stand-alone library for CHIMES Driver, but it can also be used as a guide for how to add the CHIMES source and header files to your own Makefile. 

The sections below describe the data structures that will be needed for CHIMES, and how to initialise and call the CHIMES chemistry solver. 

.. toctree::
   :maxdepth: 1

   Data_Structures
   Init_Chimes
   Call_Chimes
   Parallelisation
   Redshift_UVB
   Hybrid_Cool
   Misc
