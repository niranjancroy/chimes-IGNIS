.. CHIMES Driver Inputs
   Alexander Richings, 4th March 2020

.. _ChimesDriverInputs_label: 

Inputs
------

There are two ways in which we can define the input physical parameters of the gas that we want to apply the CHIMES chemistry solver to. We can either read in the gas particle data from an HDF5 snapshot file (using ``IO_mode == snapshot``), or we can define a grid of temperatures, densities and metallicities (using ``IO_mode == grid``). These options are described in more detail below. 

Snapshot
^^^^^^^^

You can provide an HDF5 snapshot from a simulation, as specified by the ``input_file`` parameter. CHIMES Driver will then read in the gas particles from the snapshot and integrate the chemistry for each particle. Note that these snapshot files don't actually have to be outputs from a simulation. They can also be used to define any general gas distribution by hand that can then be passed to CHIMES Driver. 

We have provided routines to handle snapshots in the format produced by GIZMO or AREPO. This can be set via the ``snapshot_type`` parameter. This will automatically find most of the data arrays such as gas densities as appropriate for the given snapshot type. For some of the arrays the user can specify the name of the array in the snapshot using certain parameters in the parameter file (see the :ref:`ChimesDriverParam_label` section for details), but in general these can be left to their default values. 

Alternatively, you can provide your own snapshot load function, if your snapshot type is not already supported. If you set ``snapshot_type = USER``, it will then look for a Python module in ``utils/user_def_utils.py``. This file is not provided in the ``chimes-driver.py`` repository; instead, the user will need to create their own copy of this file to define their own routines. 

It will then look for a routine called ``user_load_function(my_snapshot_data)`` in this Python module file. This routine must take an instance of the ``SnapshotData`` class (defined in ``utils/snapshot_utils.py``) as an argument. This class contains the various particle arrays that will need to be read in from the snapshot. 

You can look at the ``load_GIZMO(self)`` and ``load_AREPO(self)`` routines in ``utils/snapshot_utils.py`` for examples of how to read the snapshot data into the ``SnapshotData`` class. Not all of the data arrays are necessarily required, depending on the set up that you are using. The various data arrays are described below: 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``nH_arr``                          | | Hydrogen number density, in units of cm^-3. Always required.               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``temperature_arr``                 | | Temperature in units of K. Always required.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``metallicity_arr``                 | | Two-dimensional array of size (``N_gas`` x 11), where ``N_gas`` is the     |
|                                     | | number of gas particles, containing the mass fraction of each element      |
|                                     | | relative to the total mass, in the order: All metals, He, C, N, O, Ne,     |
|                                     | | Mg, Si, S, Ca, Fe. Always required.                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``shieldLength_arr``                | | Shielding length, in units of cm. This does not need to be read in by the  |
|                                     | | ``user_load_function(my_snapshot_data)`` routine, because it will          |
|                                     | | instead be set by the ``set_shielding_array(self)`` routine in             |
|                                     | | in ``utils/snapshot_utils.py``, according to the ``shield_mode``.          |
|                                     | | If you are reading in the shielding lengths from the snapshot, i.e. with   |
|                                     | | ``shield_mode == read-in``, you just need to make sure that the            |
|                                     | | ``snapshot_column_density_array`` parameter points to the correct array    |
|                                     | | in the snapshot. The shielding lengths are then calculated by dividing     |
|                                     | | these column densities by ``nH_arr``.                                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``init_chem_arr``                   | | Two-dimensional array of size (``N_gas`` x ``N_species``), where           |
|                                     | | ``N_gas`` is the number of gas particles and ``N_species`` is the number   |
|                                     | | of species in the network. This only needs to be read in from the          |
|                                     | | snapshot if ``snapshot_chemistry_array != None``.                          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``ChimesFluxIon_arr``               | | Two-dimensional array of size (``N_gas`` x 8), where ``N_gas`` is the      |
|                                     | | number of gas particles, giving the fluxes in the >13.6 eV band in the     |
|                                     | | 8 stellar age bins as used in Richings et al. (in prep). These fluxes      |
|                                     | | are given as ``isotropic_photon_density * speed_of_light`` in units of     |
|                                     | | cm^-2 s^-1. These only need to be read in from the snapshot if             |
|                                     | | ``UV_field == StellarFluxes`` and ``compute_stellar_fluxes == 0``.         |
|                                     | | In the ``load_GIZMO(self)`` and ``load_AREPO(self)`` routines, the name    |
|                                     | | of this array in the snapshot is defined by the                            |
|                                     | | ``snapshot_flux_ion_array`` parameter. If we instead compute the fluxes    |
|                                     | | in CHIMES Driver, i.e. with ``compute_stellar_fluxes == 1``, then this     |
|                                     | | array will be written out to ``snapshot_flux_ion_array`` in the output     |
|                                     | | HDF5 file.                                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``ChimesFluxG0_arr``                | | Two-dimensional array of size (``N_gas`` x 8), where ``N_gas`` is the      |
|                                     | | number of gas particles, giving the fluxes in the 6-13.6 eV band in the    |
|                                     | | 8 stellar age bins as used in Richings et al. (in prep). These fluxes      |
|                                     | | are given as the strength of the radiation field in Habing units. These    |
|                                     | | only need to be read in from the snapshot if                               |
|                                     | | ``UV_field == StellarFluxes`` and ``compute_stellar_fluxes == 0``.         |
|                                     | | In the ``load_GIZMO(self)`` and ``load_AREPO(self)`` routines, the name    |
|                                     | | of this array in the snapshot is defined by the                            |
|                                     | | ``snapshot_flux_G0_array`` parameter. If we instead compute the fluxes     |
|                                     | | in CHIMES Driver, i.e. with ``compute_stellar_fluxes == 1``, then this     |
|                                     | | array will be written out to ``snapshot_flux_G0_array`` in the output      |
|                                     | | HDF5 file.                                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``gas_coords_arr``                  | | Two-dimensional array of size (``N_gas`` x 3), where ``N_gas`` is the      |
|                                     | | number of gas particles, giving the coordinates of the gas particles in    |
|                                     | | units of cm. Used to calculate the stellar fluxes. Only needs to be read   |
|                                     | | in from the snapshot if ``UV_field == StellarFluxes`` and                  |
|                                     | | ``compute_stellar_fluxes == 1``.                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``star_coords_arr``                 | | Two-dimensional array of size (``N_star`` x 3), where ``N_star`` is the    |
|                                     | | number of star particles, giving the coordinates of the star particles in  |
|                                     | | units of cm. Used to calculate the stellar fluxes. Only needs to be read   |
|                                     | | in from the snapshot if ``UV_field == StellarFluxes`` and                  |
|                                     | | ``compute_stellar_fluxes == 1``.                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``star_mass_arr``                   | | The masses of the star particles in units of g. Used to calculate the      |
|                                     | | stellar fluxes. Only needs to be read in from the snapshot if              |
|                                     | | ``UV_field == StellarFluxes`` and ``compute_stellar_fluxes == 1``.         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``star_age_Myr_arr``                | | The stellar ages of the star particles in units of Myr. Used to calculate  |
|                                     | | the stellar fluxes. Only needs to be read in from the snapshot if          |
|                                     | | ``UV_field == StellarFluxes`` and ``compute_stellar_fluxes == 1``.         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``HIIregion_delay_time``            | | For gas particles that are flagged as HII regions, this arrays gives how   |
|                                     | | long it will remain as an HII region. This array is zero for particles     |
|                                     | | that are not in an HII region, and so can be used to identify particles    |
|                                     | | that have been flagged as an HII region. The name of this array in the     |
|                                     | | snapshot is given by the ``snapshot_HIIregion_array`` parameter. This is   |
|                                     | | only used if ``disable_shielding_in_HII_regions == 1``.                    | 
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

Grid
^^^^

You can also set up a grid of gas temperatures, densities and metallicities and the evolve the chemistry on each point in the grid. This is particularly useful for creating tables of equilibrium abundances or cooling rates, but can also be used to explore the non-equilibrium chemical evolution across a regular grid of physical conditions. The parameters that are used to define the minimum, maximum and logarithmic spacing in each of the grid dimensions are described in the :ref:`ChimesDriverParam_label` section, under the **General Parameters**. 
