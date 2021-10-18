.. CHIMES Driver Parameters 
   Alexander Richings, 3rd March 2020 

.. _ChimesDriverParam_label: 

Parameters
----------

The following parameters can be set via a text file, which is then passed to the ``chimes-driver.py`` script as an argument. If a parameter is not specified in this file, then it will use the default value, as defined in ``chimes-driver/driver_config.py``. 

General Parameters
^^^^^^^^^^^^^^^^^^

The following parameters control the general behaviour of the Driver script. 

+-------------------------------------+------------------------------------------------------------------------------+
| Parameter                           | Description                                                                  |
+=====================================+==============================================================================+
| ``chimes_library_path``             | | Path to the ``libchimes.so`` library file created when you built the       |
|                                     | | CHIMES module (see the :ref:`build_label` section).                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``chimes_data_path``                | | Path to the ``chimes-data`` repository that you cloned from Bit Bucket     |
|                                     | | (see the :ref:`Download_label` section).                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``EqAbundanceTable_``             | | Path to the equilibrium abundance table, relative to the                   |
| | ``filename``                      | | ``EqAbundancesTables`` directory in the ``chimes-data`` repository. If you |
|                                     | | are not using equilibrium cooling, you can leave this as the default       |
|                                     | | value.                                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IO_mode``                         | | Determines how the input gas properties will be defined. Possible values   |
|                                     | | are as follows:                                                            |
|                                     | | ``snapshot`` - Read in gas particle data from an HDF5 snapshot.            |
|                                     | | ``grid`` - Set up a grid of gas temperatures, densities and metallicities. |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``driver_mode``                     | | Determines the mode of operation that will be used. Possible values are    |
|                                     | | follows:                                                                   |
|                                     | | ``eqm_state`` - evolve the chemistry for a period of time and record the   |
|                                     | |   final abundances as ``n_i / n_Htot``.                                    |
|                                     | | ``eqm_table`` - as ``eqm_state``, but records the final abundances         |
|                                     | |   relative to the corresponding element, i.e. ``n_i / n_elem``. This is    |
|                                     | |   used to create the equilibrium abundance tables. Can only be used with   |
|                                     | |   ``IO_mode == grid``.                                                     |
|                                     | | ``cooling_rates`` - evolve the chemistry for a period of time and          |
|                                     | |   record the total cooling and heating rates. These are given as           |
|                                     | |   log10(rate [erg cm^-3 s^-1]).                                            |
|                                     | | ``noneq_evolution`` - evolve the chemistry for a period of time, and       |
|                                     | |   record the evolution of the abundances (as ``n_i / n_Htot``) and         |
|                                     | |   the temperature.                                                         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``UV_field``                        | | Sets the UV radiation field that will be used in the chemistry solver.     |
|                                     | | Possible options are: ``None``, ``HM01``, ``HM12``, ``FG20``, ``B87``      |
|                                     | | ``Colibre``, ``StellarFluxes`` and ``S04``. See the                        |
|                                     | | :ref:`ChimesDriverUVField_label` section for more details.                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``shield_mode``                     | | Sets the Shielding mode that will be used in the chemistry solver          |
|                                     | | Possible options are: ``None``, ``read-in``, ``Jeans``, ``Colibre``. See   |
|                                     | | the :ref:`ChimesDriverShielding_label` section for more details.           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``dust_depletion``                  | | Sets how the depletion of metals on to dust grains will be treated in      |
|                                     | | the chemistry solver. Possible values are: ``None``, ``J09``, ``DC16`` and |
|                                     | | ``Colibre``. See the :ref:`ChimesDriverDustDepletion_label` section for    |
|                                     | | more details.                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``input_file``                      | | Name of the HDF5 snapshot file to be used as input. This parameter         |
|                                     | | is only used if ``IO_mode == snapshot``.                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``output_file``                     | | Name of the HDF5 file that the outputs will be written to. If              |
|                                     | | ``IO_mode == snapshot``, you can set the output file to be the same as     |
|                                     | | the input file, in which the output data arrays will be appended to this   |
|                                     | | this file, or you can specify a different file. The script will check      |
|                                     | | that, if the output file already exists, the output data arrays are not    |
|                                     | | already present, otherwise it will abort rather than try to overwrite      |
|                                     | | them.                                                                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``hdf5_output_group``               | | Name of the HDF5 group within the output file where the output arrays      |
|                                     | | will be written out to. If this is not given in the parameter file,        |
|                                     | | it defaults to "/", i.e. the root of the HDF5 file. This is only used if   |
|                                     | | ``IO_mode == snapshot``. See the :ref:`ChimesDriverOutputs_label` section  |
|                                     | | for details.                                                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``snapshot_type``                   | | Specifies which simulation code produced the input snapshot file.  Only    |
|                                     | | used when ``IO_mode == snapshot``. Possible options are: ``None``,         |
|                                     | | ``GIZMO``, ``AREPO`` and ``USER``. See the :ref:`ChimesDriverInputs_label` |
|                                     | | section for more details.                                                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``snapshot_cosmo_flag``             | | Integer flag that determines whether the input snapshot file was a         |
|                                     | | cosmological run, i.e. with co-moving units. Only used when                |
|                                     | | ``IO_mode == snapshot`` and for certain ``snapshot_type`` values.          |
|                                     | | Possible options are ``0`` (non-cosmological) and ``1`` (cosmological).    |
|                                     | | See the :ref:`ChimesDriverInputs_label` section for more details.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``snapshot_``                     | | Unit mass in cgs units used in the snapshots. Only used when               |
| | ``unitMass_cgs``                  | | ``IO_mode == snapshot``. Defaults to ``1.989e43``, i.e. 10^10 Msol.        |
|                                     | | See the :ref:`ChimesDriverInputs_label` section for more details.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``snapshot_``                     | | Unit length in cgs units used in the snapshots. Only used when             |
| | ``unitLength_cgs``                | | ``IO_mode == snapshot``. Defaults to ``3.0857e21``, i.e. 1 kpc. See the    |
|                                     | | :ref:`ChimesDriverInputs_label` section for more details.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``snapshot_``                     | | Unit velocity in cgs units used in the snapshots. Only used when           |
| | ``unitVelocity_cgs``              | | ``IO_mode == snapshot``. Defaults to ``1.0e5``, i.e. 1 km/s. See the       |
|                                     | | :ref:`ChimesDriverInputs_label` section for more details.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``snapshot_``                     | | Name of the array in the input HDF5 snapshot file that contains the        |
| | ``chemistry_array``               | | initial CHIMES abundance array for each gas particle. If this is set to    |
|                                     | | ``None``, we will instead set the initial abundances to an initial state   |
|                                     | | as determined by the ``InitIonState`` parameter (see below). If this is    |
|                                     | | set to a value other than ``None`` but the specified array is not          |
|                                     | | present in the snapshot, it will exit with an error message. See the       |
|                                     | | :ref:`ChimesDriverInputs_label` section for more details.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``snapshot_column_``              | | Name of the array in the input HDF5 snapshot file that contains the        |
| | ``density_array``                 | | hydrogen column density of each gas particle that will be used to shield   |
|                                     | | the UV radiation field. Only used if ``IO_mode == snapshot`` and           |
|                                     | | ``shield_mode == read-in``. See the :ref:`ChimesDriverInputs_label`        |
|                                     | | section for more details.                                                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``snapshot_flux_``                | | If ``compute_stellar_fluxes == 0``, this is the name of the array in the   |
| | ``ion_array``                     | | input HDF5 snapshot file that contains the stellar fluxes in the           |
|                                     | | >13.6 eV band.                                                             |
|                                     | | If ``compute_stellar_fluxes == 1``, the stellar fluxes in the >13.6 eV     |
|                                     | | band will be written out to this array in the output HDF5 file, unless     |
|                                     | | this parameter is set to ``None``, in which case the fluxes are not        |
|                                     | | written out.                                                               |
|                                     | | Only used if ``IO_mode == snapshot`` and ``UV_field == StellarFluxes``.    |
|                                     | | See the :ref:`ChimesDriverInputs_label` and                                |
|                                     | | :ref:`ChimesDriverOutputs_label` sections for more details.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``snapshot_flux_``                | | If ``compute_stellar_fluxes == 0``, this is the name of the array in the   |
| | ``G0_array``                      | | input HDF5 snapshot file that contains the stellar fluxes in the           |
|                                     | | 6-13.6 eV band.                                                            |
|                                     | | If ``compute_stellar_fluxes == 1``, the stellar fluxes in the 6-13.6 eV    |
|                                     | | band will be written out to this array in the output HDF5 file, unless     |
|                                     | | this parameter is set to ``None``, in which case the fluxes are not        |
|                                     | | written out.                                                               |
|                                     | | Only used if ``IO_mode == snapshot`` and ``UV_field == StellarFluxes``.    |
|                                     | | See the :ref:`ChimesDriverInputs_label` and                                |
|                                     | | :ref:`ChimesDriverOutputs_label` sections for more details.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``compute_stellar_``              | | Integer flag that determines whether to compute the stellar fluxes from    |
| | ``fluxes``                        | | the star particles in the snapshot. Only used if ``IO_mode == snapshot``   |
|                                     | | and ``UV_field == StellarFluxes``. Possible values are ``0`` (read in      |
|                                     | | fluxes from the snapshot) or ``1`` (compute fluxes from the star           |
|                                     | | particles). See the :ref:`ChimesDriverUVField_label` section for more      |
|                                     | | details.                                                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``stellar_fluxes_``               | | Escape fraction of stellar radiation from star particles in the >13.6 eV   |
| | ``fEsc_ion``                      | | band. Only used if ``IO_mode == snapshot``, ``UV_field == StellarFluxes``  |
|                                     | | and ``compute_stellar_fluxes == 1``. See the                               |
|                                     | | :ref:`ChimesDriverUVField_label` section for more details.                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``stellar_fluxes_``               | | Escape fraction of stellar radiation from star particles in the            |
| | ``fEsc_G0``                       | | 6-13.6 eV band. Only used if ``IO_mode == snapshot``,                      |
|                                     | | ``UV_field == StellarFluxes`` and ``compute_stellar_fluxes == 1``. See     |
|                                     | | the :ref:`ChimesDriverUVField_label` section for more details.             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``disable_shielding_``            | | Integer flag to determine whether to disable shielding for particles       |
| | ``in_HII_regions``                | | flagged as HII regions. Possble values are ``0`` (do not disable           |
|                                     | | shielding) or ``1`` (disable shielding). Only used if                      |
|                                     | | ``IO_mode == snapshot``.                                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``snapshot_HIIregion_array``      | | Name of the array in the input HDF5 snapshot file that contains the        |
|                                     | | HII region delay time for each particle. Only used if                      |
|                                     | | ``disable_shielding_in_HII_regions`` == 1. See the                         |
|                                     | | :ref:`ChimesDriverInputs_label` section for more details.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``n_iterations``                    | | The number of iterations to integrate the chemistry for each gas particle  |
|                                     | | or grid point. For each iteration, it will be integrated for a period of   |
|                                     | | time set by the ``hydro_timestep`` parameter (see below). The chemistry    |
|                                     | | solver will then internally sub-cycle each iteration to perform the        |
|                                     | | integration itself. If you are using the ``eqm_state``, ``eqm_table`` or   |
|                                     | | ``cooling_rates`` Driver Modes, you would typically want to evolve the     |
|                                     | | chemistry to equilibrium. However, if you are including shielding, you     |
|                                     | | will need to do this with multiple iterations, because the column          |
|                                     | | densities of individual species (HI, H2 etc.) depend on the abundances,    |
|                                     | | but they do not get updated throughout the course of the integration. So   |
|                                     | | with each iteration, the column densities get updated according to the     |
|                                     | | current abundance array, and then the abundances are integrated at fixed   |
|                                     | | column densities. 10 iterations is usually sufficient to converge on an    |
|                                     | | equilibrium solution, but you will need to test this for your particular   |
|                                     | | set up.                                                                    |
|                                     | | If you are using the ``noneq_evolution`` Driver Mode, it will record the   |
|                                     | | chemical abundances and temperature after each iteration. The number of    |
|                                     | | iterations therefore controls how many time outputs will be recorded for   |
|                                     | | each gas particle or grid point.                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_T_min``                       | | The log10 of the minimum temperature when constructing the grid of         |
|                                     | | temperatures, densities and metallicities. Only used if                    |
|                                     | | ``IO_mode == grid``. See the :ref:`ChimesDriverInputs_label` section for   |
|                                     | | more details.                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_T_max``                       | | The log10 of the maximum temperature when constructing the grid of         |
|                                     | | temperatures, densities and metallicities. Only used if                    |
|                                     | | ``IO_mode == grid``. See the :ref:`ChimesDriverInputs_label` section for   |
|                                     | | more details.                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``delta_log_T``                     | | The logarithmic spacing between temperatures when constructing the grid    |
|                                     | | of temperatures, densities and metallicities. Only used if                 |
|                                     | | ``IO_mode == grid``. See the :ref:`ChimesDriverInputs_label` section for   |
|                                     | | more details.                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_nH_min``                      | | The log10 of the minimum hydrogen number density (in units of cm^-3)       |
|                                     | |  when constructing the grid of temperatures, densities and                 |
|                                     | | metallicities. Only used if ``IO_mode == grid``. See the                   |
|                                     | | :ref:`ChimesDriverInputs_label` section for more details.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_nH_max``                      | | The log10 of the maximum hydrogen number density (in units of cm^-3)       |
|                                     | | when constructing the grid of temperatures, densities and                  |
|                                     | | metallicities. Only used if ``IO_mode == grid``. See the                   |
|                                     | | :ref:`ChimesDriverInputs_label` section for more details.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``delta_log_nH``                    | | The logarithmic spacing between densities when constructing the grid of    |
|                                     | | temperatures, densities and metallicities. Only used if                    |
|                                     | | ``IO_mode == grid``. See the :ref:`ChimesDriverInputs_label` section for   |
|                                     | | more details.                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_Z_min``                       | | The log10 of the minimum metallicity (relative to solar metallicity,       |
|                                     | | Zsol = 0.0129) when constructing the grid of temperatures, densities and   |
|                                     | | metallicities. Only used if ``IO_mode == grid``. See the                   |
|                                     | | :ref:`ChimesDriverInputs_label` section for more details.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_Z_max``                       | | The log10 of the maximum metallicity (relative to solar metallicity,       |
|                                     | | Zsol = 0.0129) when constructing the grid of temperatures, densities and   |
|                                     | | metallicities. Only used if ``IO_mode == grid``. See the                   |
|                                     | | :ref:`ChimesDriverInputs_label` section for more details.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``delta_log_Z``                     | | The logarithmic spacing between metallicities when constructing the grid   |
|                                     | | of temperatures, densities and metallicities. Only used if                 |
|                                     | | ``IO_mode == grid``. See the :ref:`ChimesDriverInputs_label` section for   |
|                                     | | more details.                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``shield_length_``                | | Factor to multiply the shielding length by for all gas particles/grid      |
| | ``factor``                        | | points. Only used if ``shield_mode`` is set to ``Jeans`` or ``Colibre``,   |
|                                     | | or if either ``UV_field`` or ``dust_depletion`` are set to ``Colibre``     |
|                                     | | (because the reference column density, ``N_ref``, that is used to scale    |
|                                     | | the interstellar radiation field and the dust depletion factors in the     |
|                                     | | COLIBRE model is also multiplied by this factor; See Ploeckinger & Schaye  |
|                                     | | (2020). See the :ref:`ChimesDriverShielding_label` section for more        |
|                                     | | details.                                                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``max_shield_length``               | | Maximum shielding length, in cm. If not specified in the parameter file,   |
|                                     | | it defaults to ``3.086e23``, i.e. 100 kpc. Only used if ``shield_mode``    |
|                                     | | is set to ``Jeans`` or ``Colibre``, or if either ``UV_field`` or           |
|                                     | | ``dust_depletion`` are set to ``Colibre`` (because the reference column    |
|                                     | | density, ``N_ref``, that is used to scale the interstellar radiation       |
|                                     | | field and the dust depletion factors in the COLIBRE model is also          |
|                                     | | limited by this maximum; See Ploeckinger & Schaye 2020). See the           |
|                                     | | :ref:`ChimesDriverShielding_label` section for more details.               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``colibre_log_T_min``               | | The log_T_min parameter in the temperature scaling of the ``N_ref``        |
|                                     | | column density used in the COLIBRE model. See Ploeckinger & Schaye         |
|                                     | | (2020) for details. Only used if ``UV_field``, ``shield_mode`` or          |
|                                     | | ``dust_depletion`` are set to ``Colibre``.                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``colibre_log_T_max``               | | The log_T_max parameter in the temperature scaling of the ``N_ref``        |
|                                     | | column density used in the COLIBRE model. See Ploeckinger & Schaye         |
|                                     | | (2020) for details. Only used if ``UV_field``, ``shield_mode`` or          |
|                                     | | ``dust_depletion`` are set to ``Colibre``.                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``colibre_scale_``                | | The strength of the Milky Way interstellar radiation field used in the     |
| | ``MW_ISRF``                       | | Colibre model is scaled by this value (see Ploeckinger & Schaye 2020       |
|                                     | | for details).                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``radiation_field_``              | | The strength of the radiation field from the cross sections tables are     |
| | ``normalisation_``                | | multiplied by this factor. Only used if ``UV_field`` is set to ``HM01``,   |
| | ``factor``                        | | ``HM12``, ``FG20`` or ``B87``.                                             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``bolometric_AGN_``               | | The bolometric luminosity of the AGN, in units of erg/s. Only used if      |
| | ``luminosity_cgs``                | | ``UV_field`` is set to ``S04``.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``distance_to_AGN_kpc``           | | Distance to the AGN in kpc, which determines the strength of the           |
|                                     | | quasar UV spectrum. Only used if ``UV_field`` is set to ``S04``            |
|                                     | | and ``IO_mode`` is set to ``grid``.                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``AGN_position_x_kpc``            | | The x position of the AGN in kpc, used to determine the distance from      |
|                                     | | each gas particle in the snapshot to the AGN. Only used if ``UV_field``    |
|                                     | | is set to ``S04`` and ``IO_mode`` is set to ``snapshot``.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``AGN_position_y_kpc``            | | The y position of the AGN in kpc, used to determine the distance from      |
|                                     | | each gas particle in the snapshot to the AGN. Only used if ``UV_field``    |
|                                     | | is set to ``S04`` and ``IO_mode`` is set to ``snapshot``.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``AGN_position_z_kpc``            | | The z position of the AGN in kpc, used to determine the distance from      |
|                                     | | each gas particle in the snapshot to the AGN. Only used if ``UV_field``    |
|                                     | | is set to ``S04`` and ``IO_mode`` is set to ``snapshot``.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

Global Variables
^^^^^^^^^^^^^^^^

The following parameters control the general behaviour of the chemistry solver. They correspond to the variables stored in the **globalVariables** data structure within the CHIMES module. Note that not all of the variables in this data structure can be set directly from the parameter file in this way. 

+-------------------------------------+------------------------------------------------------------------------------+
| Parameter                           | Description                                                                  |
+=====================================+==============================================================================+
| ``redshift``                        | | Redshift used for the redshift-dependent extragalactic UV background.      |
|                                     | | Only used if ``UV_field`` is set to ``HM12``, ``FG20`` or ``Colibre``.     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reionisation_redshift``           | | Redshift of reionisation, used for the redshift-dependent extragalactic    |
|                                     | | UV background, which is disabled if the redshift is before reionisation.   |
|                                     | | Only used if ``UV_field`` is set to ``HM12``, ``FG20`` or ``Colibre``.     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``use_redshift_``                 | | Integer flag (``0`` or ``1``) that determines whether to use redshift-     |
| | ``dependent_eqm_tables``          | | dependent equilibrium abundance tables. Only used if ``UV_field`` is set   |
|                                     | | to ``HM12``, ``FG20`` or ``Colibre``.                                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``StaticMolCooling``                | | Integer flag that controls how the effective CO and H2O column densities   |
|                                     | | are calculated for the CO and H2O molecular cooling rates. Possible        |
|                                     | | options are:                                                               |
|                                     | | ``0`` - Include the divergence of the velocity in the effective column     |
|                                     | |   densities (see equation 3.1 in Richings et al 2014a). If you use this    |
|                                     | |   option, you will need to set the velocity divergence, ``divVel``, in     |
|                                     | |   the Gas Variables (see below).                                           |
|                                     | | ``1`` - Assume a static gas distribution, i.e. calculate the line          |
|                                     | |   broadening based only on the thermal temperature.                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``T_mol``                           | | Temperature cut above which the molecular network is disabled and the      |
|                                     | | molecule abundances are set to zero.                                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``grain_temperature``               | | Temperature of the dust grains, used in the rate of H2 formation on dust   |
|                                     | | and in the cooling rate due to energy transfer between gas and dust.       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``cmb_temperature``                 | | Temperature of the Cosmic Microwave Background, used for Compton           |
|                                     | | cooling from the CMB. Note that this is not automatically set according    |
|                                     | | to the specified redshift. The User will need to set this to the           |
|                                     | | required value.                                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``relativeTolerance``               | | Relative tolerance used to control the accuracy of the chemistry and       |
|                                     | | cooling integration.                                                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``absoluteTolerance``               | | Absolute tolerance used to control the accuracy of the chemistry           |
|                                     | | integration.                                                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``explicitTolerance``               | | Each time we call the CHIMES chemistry solver routine, we first calculate  |
|                                     | | the final chemical state and temperature using an explicit integration     |
|                                     | | scheme. If the relative change in the temperature and all species above    |
|                                     | | the ``absoluteTolerance`` is below the ``explicitTolerance``, we just      |
|                                     | | take this explicit solution. Otherwise, we use the full CVODE implicit     |
|                                     | | integration scheme. This allows us to avoid a lot of the expensive         |
|                                     | | overheads involved in the implicit solver when this is not needed, for     |
|                                     | | example if the particle has already reach equilibrium.                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``scale_metal_``                  | | Integer flag that controls how we set the absolute tolerances for metals.  |
| | ``tolerances``                    | | Possible values are as follows:                                            |
|                                     | | ``0`` - Use a constant absolute tolerance for all species, given by the    |
|                                     | |   ``absoluteTolerance`` parameter.                                         |
|                                     | | ``1`` - Scale the ``absoluteTolerance`` parameter for each species by      |
|                                     | |   the total abundance of its corresponding element (including He). This is |
|                                     | |   particularly important for cosmological simulations where gas particles  |
|                                     | |   can have very small element abundances. If the total abundance of a      |
|                                     | |   metal if less than ``absoluteTolerance``, then all of its ions and       |
|                                     | |   molecules will be below this tolerance and so will carry little or no    |
|                                     | |   weight in the error estimation when determining how to sub-cycle the     |
|                                     | |   integration. This can make the integration unstable, as then the sum of  |
|                                     | |   the ions and molecules of that metal is not well constrained, which can  |
|                                     | |   lead to negative abundances. By scaling the metal tolerances in this     |
|                                     | |   way we avoid this problem.                                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``chimes_debug``                    | | Integer flag that controls how much debug information is printed when we   |
|                                     | | encounter a CVODE error or warning flag. The possible values are as        |
|                                     | | follows:                                                                   |
|                                     | | ``0`` - Ignore all CVODE error and warning messages.                       |
|                                     | | ``1`` - Print only the CVODE error or warning message.                     |
|                                     | | ``2`` - Print the CVODE error or warning message, and also all of the      |
|                                     | |         variables in the **gasVariables** structure, so that we can see    |
|                                     | |         the gas properties where this error occurred.                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``hybrid_cooling_mode``             | | Integer flag (``0`` or ``1``) that controls whether to use the hybrid      |
|                                     | | cooling mode in CHIMES, where a user-defined function is used to           |
|                                     | | calculate an additional cooling rate that is added on to the cooling and   |
|                                     | | heating rates calculated in CHIMES. This can be used, for example, to add  |
|                                     | | on the cooling rates from elements that have been switched off in the      |
|                                     | | non-equilibrium network (see below) using pre-computed cooling tables.     |
|                                     | | Note that this option has not yet been implemented in                      |
|                                     | | ``chimes-driver.py``, so this parameter should always be set to ``0``.     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeCarbon``                   | | Integer flag (``0`` or ``1``) that controls whether to include carbon      |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeNitrogen``                 | | Integer flag (``0`` or ``1``) that controls whether to include nitrogen    |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeOxygen``                   | | Integer flag (``0`` or ``1``) that controls whether to include oxygen      |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeNeon``                     | | Integer flag (``0`` or ``1``) that controls whether to include neon        |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeMagnesium``                | | Integer flag (``0`` or ``1``) that controls whether to include magnesium   |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeSilicon``                  | | Integer flag (``0`` or ``1``) that controls whether to include silicon     |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeSulphur``                  | | Integer flag (``0`` or ``1``) that controls whether to include sulphur     |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeCalcium``                  | | Integer flag (``0`` or ``1``) that controls whether to include calcium     |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``IncludeIron``                     | | Integer flag (``0`` or ``1``) that controls whether to include iron        |
|                                     | | in the CHIMES network.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

Gas Variables
^^^^^^^^^^^^^

The following parameters control more specific behaviour of the chemistry solver. They correspond to the variables stored in the **gasVariables** data structure within the CHIMES module. Note that not all of the variables in this data structure can be set directly from the parameter file in this way. 

+-------------------------------------+------------------------------------------------------------------------------+
| Parameter                           | Description                                                                  |
+=====================================+==============================================================================+
| ``cr_rate``                         | | Ionisation rate of HI due to cosmic rays. Note that, if                    |
|                                     | | ``UV_field == Colibre``, this rate will then be multiplied by              |
|                                     | | ``colibre_scale_MW_ISRF * ((N_ref / N_H0) ** 1.4)`` (see Ploeckinger &     |
|                                     | | Schaye 2020).                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TempFloor``                       | | Temperature floor used in the cooling integration, in K.                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``divVel``                          | | Divergence of the velocity field, in cgs units. Only used if               |
|                                     | | ``StaticMolCooling == 0``.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``doppler_broad``                   | | Doppler broadening parameter due to turbulence, in km/s, used in the H2    |
|                                     | | self-shielding function (see Richings et al. 2014b).                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``ForceEqOn``                       | | Integer flag (``0`` or ``1``) that controls whether to set the             |
|                                     | | chemical abundances to equilibrium from the pre-computed equilibrium       |
|                                     | | abundance tables. If ``ThermEvolOn == 1``, the cooling will be evolved     |
|                                     | | using cooling and heating rates in chemical equilibrium, where it will     |
|                                     | | set the abundances to equilibrium from the tables and then compute the     |
|                                     | | cooling and heating rates from those abundances.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``ThermEvolOn``                     | | Integer flag (``0`` or ``1``) that controls whether to evolve the          |
|                                     | | temperature along with the chemistry. If set to ``0``, the chemical        |
|                                     | | abundances will be evolved at fixed temperature.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``temp_floor_mode``                 | | Integer flag that controls how we implement the temperature floor. Only    |
|                                     | | used if ``ThermEvolOn == 1``. The possible values are as follows:          |
|                                     | | ``0`` - If the temperature is less than ``TempFloor`` the rate of change   |
|                                     | |         of internal energy, ``du_dt``, is constrained to be positive,      |
|                                     | |         so that it can heat up again but it cannot cool any further. The   |
|                                     | |         chemical abundances continue to be integrated as normal.           |
|                                     | | ``1`` - If the temperature is less than ``TempFloor``, we immediately      |
|                                     | |         halt the CVODE integration. The chemical abundances and            |
|                                     | |         temperature are then kept at the values that they reached at the   |
|                                     | |         point where the integration was halted. This tends to be more      |
|                                     | |         stable. However, it means that the chemical abundances do not      |
|                                     | |         continue to evolve at a fixed temperature at the ``TempFloor``.    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``InitIonState``                    | | Integer that controls the initial chemical state that the chemistry        |
|                                     | | integration will start from. Each element will be set to the ionisation    |
|                                     | | state given by this parameter. For example, if set to ``0`` all elements   |
|                                     | | will be neutral, whereas if set to ``1`` all elements will be singly       |
|                                     | | ionised. If ``IO_mode == snapshot`` and                                    |
|                                     | | ``snapshot_chemistry_array != None``, we will instead read in the inital   |
|                                     | | chemical state from the snapshot.                                          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``constant_heating_``             | | Constant heating rate (positive for heating, negative for cooling), in     |
| | ``rate``                          | | units of erg cm^-3 s^-1, that is added on to the heating and cooling       |
|                                     | | rates from the CHIMES network.                                             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``hydro_timestep``                  | | The length of time to integrate the chemistry and cooling for each         |
|                                     | | iteration, in seconds.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

References
^^^^^^^^^^

| `Ploeckinger & Schaye (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200614322P/abstract>`_
| `Richings et al. (2014a) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.440.3349R>`_
| `Richings et al. (2014b) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.2780R>`_
