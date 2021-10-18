.. CHIMES Driver Examples
   Alexander Richings, 12th March 2020

.. _ChimesDriverExamples_label:

Examples
--------

The ``examples/`` directory in the ``chimes-driver`` repository contains several example parameter files for Chimes Driver. We describe the parameters used in each example below. 

cooling_rates
^^^^^^^^^^^^^

Creates a table of cooling and heating rates on a grid of temperature, density and metallicity. 

.. code-block:: none 

  chimes_library_path        /home/username/chimes/libchimes.so 
  chimes_data_path           /home/username/chimes-data 

Defines the paths to the CHIMES library and CHIMES data repository. 

.. code-block:: none 

  driver_mode  cooling_rates
  IO_mode      grid 

Defines the mode in which Chimes Driver will be run, i.e. to compute cooling and heating rates from the final chemical state on a grid of temperature, density and metallicity. 

.. code-block:: none 

  output_file  ./cooling_rates.hdf5 

The name of the file that the outputs will be written out to. 

.. code-block:: none 

  UV_field        B87 
  shield_mode     Jeans 
  dust_depletion  None

Sets the UV radiation field, shielding length and dust depletion model. In this example, we use the local interstellar radiation field from Black (1987), a shielding length equal to the Jeans length, and a fixed dust-to-metals ratio without accounting for the reduction of metal abundances from the gas phase due to depletion on to dust grains.

.. code-block:: none 

  log_T_min     2.0  
  log_T_max     9.0
  delta_log_T   0.1
  log_nH_min    -4.0 
  log_nH_max    4.0 
  delta_log_nH  0.5
  log_Z_min     -2.0 
  log_Z_max     0.0
  delta_log_Z   0.5 

Defines the extent and logarithmic spacing of the grid in temperature, density and metallicity. 

.. code-block:: none 

  shield_length_factor  1.0 
  max_shield_length     3.086e23 

Defines a multiplicative factor to multiply all shielding lengths by (in this case unity, i.e. we will just be using the Jeans length), and a maximum shielding length of 100 kpc. 

.. code-block:: none 

  n_iterations  10 

Defines how many iterations to integrate the chemistry for on each point in the grid. The column densities of individual species are updated between each iteration. 
 
.. code-block:: none 

  StaticMolCooling  1 

The CO and H2O cooling rates are calculated for a static gas distribution, i.e. based only on thermal line broadening. 

.. code-block:: none 

  T_mol  1.0e5  

Maximum temperature for the molecular network. Above this temperature, molecule abundances are set to zero and molecular reactions are excluded from the network. 

.. code-block:: none 

  InitIonState  1 

The chemistry integration will be started from a singly ionised initial state. 

.. code-block:: none 

  grain_temperature  10.0

The temperature of dust grains is set to 10 K. 

.. code-block:: none 

  cmb_temperature  2.725

The temperature of the Cosmic Microwave Background is set to 2.725 K. 

.. code-block:: none 

  relativeTolerance         1.0e-4 
  absoluteTolerance         1.0e-10 
  explicitTolerance         0.05 
  scale_metal_tolerances    1 

Defines the various tolerance parameters used to control the accuracy of the chemistry integration. 

.. code-block:: none 

  IncludeCarbon     1 
  IncludeNitrogen   1 
  IncludeOxygen     1 
  IncludeNeon       1 
  IncludeMagnesium  1
  IncludeSilicon    1
  IncludeSulphur    1
  IncludeCalcium    1
  IncludeIron       1

Defines which metal elements to include in the non-equilibrium network. In this example, all metals are included. 

.. code-block:: none 

  cr_rate  1.8e-16 

The ionisation rate of HI due to cosmic rays, in units of s^-1. 

.. code-block:: none 

  doppler_broad  7.1

Doppler broadening parameter for the H2 self-shielding function, in units of km/s. 

.. code-block:: none 

  ForceEqOn  0

Setting this to zero tells CHIMES to integrate the chemistry in non-equilibrium. 

.. code-block:: none 

  ThermEvolOn  0

Setting this to zero disables the temperature evolution, i.e. the temperature is held fixed as the chemical abundances are evolved. 

.. code-block:: none 

  hydro_timestep  3.16e16

The total time to integrate the chemistry over for each iteration, in seconds. This value correspond to 1 Gyr. 

eqm_table
^^^^^^^^^

Creates a table of equilibrium abundances on a grid of temperature, density and metallicity. The output HDF5 file is in a format that can be used as the equilibrium abundance table for other runs with CHIMES, which will then set the abundances to equilibrium by interpolating this table if the ``ForceEqOn`` parameter is set to 1. 

.. code-block:: none 

  chimes_library_path         /home/username/chimes/libchimes.so 
  chimes_data_path            /home/username/chimes-data 

Defines the paths to the CHIMES library and CHIMES data repository. 

.. code-block:: none 

  driver_mode    eqm_table
  IO_mode        grid 

Defines the mode in which Chimes Driver will be run, i.e. to create an equilibrium abundance table on a grid of temperature, density and metallicity. 

.. code-block:: none 

  output_file  ./eqm_table.hdf5 

The name of the file that the outputs will be written out to. 

.. code-block:: none 

  UV_field        Colibre
  shield_mode     Colibre
  dust_depletion  Colibre

Sets the UV radiation field, shielding length and dust depletion model. In this example, we use the Colibre model for all three (see Ploeckinger & Schaye 2020 for details).

.. code-block:: none 

  log_T_min     2.0  
  log_T_max     9.0
  delta_log_T   0.1
  log_nH_min    -4.0 
  log_nH_max    4.0 
  delta_log_nH  0.5 
  log_Z_min     -2.0 
  log_Z_max     0.0 
  delta_log_Z   0.5  

Defines the extent and logarithmic spacing of the grid in temperature, density and metallicity. 

.. code-block:: none 

  shield_length_factor  0.5 
  max_shield_length  3.086e23 

Defines a multiplicative factor to multiply all shielding lengths by (in this case 0.5) and a maximum shielding length of 100 kpc, as used in Ploeckinger & Schaye (2020). 

.. code-block:: none 

  colibre_log_T_min      3.0 
  colibre_log_T_max      5.0 
  colibre_scale_MW_ISRF  0.1 

Parameters used in the Colibre model for the radiation field, shielding length and dust depletion (see Ploeckinger & Schaye 2020). 

.. code-block:: none 

  n_iterations  10 

Defines how many iterations to integrate the chemistry for on each point in the grid. The column densities of individual species are updated between each iteration. 

.. code-block:: none 

  StaticMolCooling  1 

The CO and H2O cooling rates are calculated for a static gas distribution, i.e. based only on thermal line broadening. 

.. code-block:: none 

  T_mol  1.0e5  

Maximum temperature for the molecular network. Above this temperature, molecule abundances are set to zero and molecular reactions are excluded from the network. 

.. code-block:: none 

  InitIonState  1 

The chemistry integration will be started from a singly ionised initial state. 

.. code-block:: none 

  grain_temperature  10.0

The temperature of dust grains is set to 10 K. 

.. code-block:: none 

  cmb_temperature  2.725

The temperature of the Cosmic Microwave Background is set to 2.725 K. 

.. code-block:: none 

  relativeTolerance         1.0e-4 
  absoluteTolerance         1.0e-10 
  explicitTolerance         0.05 
  scale_metal_tolerances    1 

Defines the various tolerance parameters used to control the accuracy of the chemistry integration. 

.. code-block:: none 

  IncludeCarbon     1 
  IncludeNitrogen   1 
  IncludeOxygen     1 
  IncludeNeon       1 
  IncludeMagnesium  1
  IncludeSilicon    1
  IncludeSulphur    1
  IncludeCalcium    1
  IncludeIron       1

Defines which metal elements to include in the non-equilibrium network. In this example, all metals are included. 

.. code-block:: none 

  redshift               0.0
  reionisation_redshift  7.5

Defines the redshift and reionisation redshift used for the extragalactic UV background component of the Colibre radiation field model. 
 
.. code-block:: none 

  cr_rate  1.8e-16 

The ionisation rate of HI due to cosmic rays, in units of s^-1. With the Colibre UV_field model, the ``cr_rate`` then gets scaled along with the interstellar radiation field component of the radiation field (see Ploeckinger & Schaye 2020). 

.. code-block:: none 

  doppler_broad  7.1

Doppler broadening parameter for the H2 self-shielding function, in units of km/s. 

.. code-block:: none 

  ForceEqOn  0

Setting this to zero tells CHIMES to integrate the chemistry in non-equilibrium. 

.. code-block:: none 

  ThermEvolOn  0

Setting this to zero disables the temperature evolution, i.e. the temperature is held fixed as the chemical abundances are evolved. 

.. code-block:: none 

  hydro_timestep  3.16e16

The total time to integrate the chemistry over for each iteration, in seconds. This value correspond to 1 Gyr. 

GIZMO_snapshot_eqm_state
^^^^^^^^^^^^^^^^^^^^^^^^

Computes the equilibrium chemical state for each gas particle in a GIZMO snapshot. 

.. code-block:: none 

  chimes_library_path         /home/username/chimes/libchimes.so 
  chimes_data_path            /home/username/chimes-data 

Defines the paths to the CHIMES library and CHIMES data repository. 

.. code-block:: none 

  driver_mode    eqm_state
  IO_mode        snapshot 

Defines the mode in which Chimes Driver will be run, i.e. to compute the equilibrium chemical abundances for each gas particle in a snapshot. 

.. code-block:: none 

  input_file   snapshot_500.hdf5
  output_file  snap500_eqm.hdf5

The name of the input snapshot file and the output file. In this example, we write the final equilibrium abundances to a separate file. 

.. code-block:: none 

  snapshot_type  GIZMO

Defines the type of snapshot. In this example, we use a GIZMO snapshot. 

.. code-block:: none 

  snapshot_chemistry_array  None

By setting this to ``None``, Chimes Driver will not try to read in the initial chemistry abundances from the snapshot. Instead, we will simply set the initial chemical state by hand, according to the ``InitIonState`` parameter (see below). 

.. code-block:: none 

  snapshot_flux_ion_array  ChimesFluxIon_output
  snapshot_flux_G0_array   ChimesFluxG0_output

Defines the names of the arrays containing the stellar fluxes for each gas particle, in the >13.6 eV and the 6-13.6 eV bands, respectively. Since we will be computing the stellar fluxes from the star particles in this example (see below), the resulting stellar fluxes will be written out to these arrays in the output file. 

.. code-block:: none 

  snapshot_cosmo_flag  0

This tells Chimes Driver that the snapshot is from a non-cosmological simulation, i.e. the units are physical, not co-moving. 

.. code-block:: none 

  snapshot_unitMass_cgs      1.989e43 
  snapshot_unitLength_cgs    3.0857e21 
  snapshot_unitVelocity_cgs  1.0e5 

Defines conversions between code units and cgs units in the snapshot. 

.. code-block:: none 

  compute_stellar_fluxes  1 

By setting this to 1, we will compute the stellar fluxes for each gas particle from the star particles in the snapshot, according to Richings et al. (in prep). 

.. code-block:: none 

  stellar_fluxes_fEsc_ion  0.05 
  stellar_fluxes_fEsc_G0   0.1 

Defines the escape fractions of radiation in each band from the neighbourhood of each star particle. 

.. code-block:: none 

  UV_field        StellarFluxes 
  shield_mode     Jeans 
  dust_depletion  DC16

Sets the UV radiation field, shielding length and dust depletion model. In this example, we use the stellar fluxes from star particles (see Richings et al. in prep), a shielding length equal to the Jeans length, and dust depletion factors based on the observations of De Cia et al. (2016). 

.. code-block:: none 

  n_iterations  10 

Defines how many iterations to integrate the chemistry for on each point in the grid. The column densities of individual species are updated between each iteration. 

.. code-block:: none 

  StaticMolCooling  1 

The CO and H2O cooling rates are calculated for a static gas distribution, i.e. based only on thermal line broadening. 

.. code-block:: none 

  T_mol  1.0e5  

Maximum temperature for the molecular network. Above this temperature, molecule abundances are set to zero and molecular reactions are excluded from the network. 

.. code-block:: none 

  InitIonState  1 

The chemistry integration will be started from a singly ionised initial state. 

.. code-block:: none 

  grain_temperature  10.0

The temperature of dust grains is set to 10 K. 

.. code-block:: none 

  cmb_temperature  2.725

The temperature of the Cosmic Microwave Background is set to 2.725 K. 

.. code-block:: none 

  relativeTolerance         1.0e-4 
  absoluteTolerance         1.0e-10 
  explicitTolerance         0.05 
  scale_metal_tolerances    1 

Defines the various tolerance parameters used to control the accuracy of the chemistry integration. 

.. code-block:: none 

  IncludeCarbon     1 
  IncludeNitrogen   1 
  IncludeOxygen     1 
  IncludeNeon       1 
  IncludeMagnesium  1
  IncludeSilicon    1
  IncludeSulphur    1
  IncludeCalcium    1
  IncludeIron       1

Defines which metal elements to include in the non-equilibrium network. In this example, all metals are included. 

.. code-block:: none 

  cr_rate  1.8e-16 

The ionisation rate of HI due to cosmic rays, in units of s^-1. 

.. code-block:: none 

  doppler_broad  7.1

Doppler broadening parameter for the H2 self-shielding function, in units of km/s. 

.. code-block:: none 

  ForceEqOn  0

Setting this to zero tells CHIMES to integrate the chemistry in non-equilibrium. 

.. code-block:: none 

  ThermEvolOn  0

Setting this to zero disables the temperature evolution, i.e. the temperature is held fixed as the chemical abundances are evolved. 

.. code-block:: none 

  hydro_timestep  3.16e16

The total time to integrate the chemistry over for each iteration, in seconds. This value correspond to 1 Gyr. 

grid_noneq_evolution
^^^^^^^^^^^^^^^^^^^^

Records the non-equilibrium evolution of temperature and chemical abundances on a grid of initial temperature, density and metallicity. 

.. code-block:: none 

  chimes_library_path        /home/username/chimes/libchimes.so 
  chimes_data_path           /home/username/chimes-data 

Defines the paths to the CHIMES library and CHIMES data repository. 

.. code-block:: none 

  driver_mode  noneq_evolution
  IO_mode      grid 

Defines the mode in which Chimes Driver will be run, i.e. to compute the non-equilibrium evolution of temperature and chemical abundances on a grid of initial temperature, density and metallicity. 

.. code-block:: none 

  output_file  ./grid_noneq_evolution.hdf5 

The name of the file that the outputs will be written out to. 

.. code-block:: none 

  UV_field        HM12 
  shield_mode     Jeans 
  dust_depletion  None

Sets the UV radiation field, shielding length and dust depletion model. In this example, we use the extragalactic UV background from Haardt & Madau (2012), a shielding length equal to the Jeans length, and a fixed dust-to-metals ratio without accounting for the reduction of metal abundances from the gas phase due to depletion on to dust grains.

.. code-block:: none 

  log_T_min     2.0  
  log_T_max     8.0
  delta_log_T   1.0
  log_nH_min    -4.0 
  log_nH_max    4.0 
  delta_log_nH  1.0
  log_Z_min     -2.0 
  log_Z_max     0.0
  delta_log_Z   1.0 

.. code-block:: none 

  shield_length_factor  1.0 
  max_shield_length     3.086e23 

Defines a multiplicative factor to multiply all shielding lengths by (in this case unity, i.e. we will just be using the Jeans length), and a maximum shielding length of 100 kpc. 

.. code-block:: none 

  n_iterations  1000 

Defines how many iterations to integrate the chemistry for on each point in the grid. We record the temperature and chemical state after each iteration. 

.. code-block:: none 

  StaticMolCooling  1 

The CO and H2O cooling rates are calculated for a static gas distribution, i.e. based only on thermal line broadening. 

.. code-block:: none 

  T_mol  1.0e5  

Maximum temperature for the molecular network. Above this temperature, molecule abundances are set to zero and molecular reactions are excluded from the network. 

.. code-block:: none 

  InitIonState  1 

The chemistry integration will be started from a singly ionised initial state. 

.. code-block:: none 

  grain_temperature  10.0

The temperature of dust grains is set to 10 K. 

.. code-block:: none 

  cmb_temperature  2.725

The temperature of the Cosmic Microwave Background is set to 2.725 K. 

.. code-block:: none 

  relativeTolerance         1.0e-4 
  absoluteTolerance         1.0e-10 
  explicitTolerance         0.05 
  scale_metal_tolerances    1 

Defines the various tolerance parameters used to control the accuracy of the chemistry integration. 

.. code-block:: none 

  IncludeCarbon     1 
  IncludeNitrogen   1 
  IncludeOxygen     1 
  IncludeNeon       1 
  IncludeMagnesium  1
  IncludeSilicon    1
  IncludeSulphur    1
  IncludeCalcium    1
  IncludeIron       1

Defines which metal elements to include in the non-equilibrium network. In this example, all metals are included. 

.. code-block:: none 

  redshift               0.0
  reionisation_redshift  7.5

Defines the redshift and reionisation redshift used for the Haardt & Madau (2012) extragalactic UV background. 

.. code-block:: none 

  cr_rate  1.8e-16 

The ionisation rate of HI due to cosmic rays, in units of s^-1. 

.. code-block:: none 

  doppler_broad  7.1

Doppler broadening parameter for the H2 self-shielding function, in units of km/s. 

.. code-block:: none 

  ForceEqOn  0

Setting this to zero tells CHIMES to integrate the chemistry in non-equilibrium. 

.. code-block:: none 

  ThermEvolOn  1

Setting this to one enables the temperature evolution. 

.. code-block:: none 

  hydro_timestep  3.16e13

The total time to integrate the chemistry over for each iteration, in seconds. This value correspond to 1 Myr. 

References
^^^^^^^^^^

| `Black (1987) <https://ui.adsabs.harvard.edu/abs/1987ASSL..134..731B>`_ 
| `De Cia et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016A%26A...596A..97D>`_
| `Haardt & Madau (2012) <https://ui.adsabs.harvard.edu/abs/2012ApJ...746..125H>`_
| `Ploeckinger & Schaye (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200614322P/abstract>`_








