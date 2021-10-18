.. Data Structures
   Alexander Richings, 18th March 2020

.. _DataStructures_label:

Data Structures
---------------

The CHIMES variables that need to be defined and set from the hydrodynamics code are broadly arranged into two structures, as described below. These are defined in the ``chimes_proto.h`` header file. 

Global Variables
^^^^^^^^^^^^^^^^

The ``globalVariables`` structure contains variables that control the overall behaviour of the CHIMES chemistry solver. You will need to create an instance of this structure within the hydrodynamics code. Most of these variables would then be set by the User, for example via a parameter file, and are then kept fixed throughout the simulation (but there are some exceptions). The variables in this structure are described below. 

+-------------------------------------+------------------------------------------------------------------------------+
| Variable                            | Description                                                                  |
+=====================================+==============================================================================+
| | ``MainDataTable``                 | | Full path to the ``chimes_main_data.hdf5`` data file (see the              |
| | ``Path[500]``                     | | :ref:`ChimesData_label` section). Note that this string can have at most   |
|                                     | | 500 characters.                                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``PhotoIonTable``                 | | Array of strings containing the full paths to each of the cross sections   |
| | ``Path[CHIMES_MAX_UV_``           | | tables, one per UV spectrum (see the :ref:`ChimesData_label` section). The |
| | ``SPECTRA][500]``                 | | maximum number of UV spectra is set to 20 in ``chimes_proto.h``. Each      |
|                                     | | path can contain at most 500 characters.                                   |
|                                     | | If you include a redshift-dependent UV background, then the path that you  |
|                                     | | give here for that spectrum needs to point to the directory containing     |
|                                     | | all of the cross sections files for each redshift, rather than the file    |
|                                     | | itself (see the :ref:`RedshiftUVB_label` section).                         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``EqAbundanceTable``              | | Full path to the equilibrium abundance table. This table is used to set    |
| | ``Path[500]``                     | | the abundance array to equilibrium when the ``ForceEqOn`` flag is set to   |
|                                     | | 1 (see the **Gas Variables** below). This string can have at most 500      |
|                                     | | characters. If you are not using the equilibrium tables at all, you can    |
|                                     | | set this string to ``None`` or ``none`` and it will skip loading the       |
|                                     | | equilibrium tables altogether.                                             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``cellSelfShieldingOn``             | | Integer flag to switch on self-shielding for each gas particle/cell.       |
|                                     | | Possible options are as follows:                                           |
|                                     | | ``0`` - Disables self-shielding. Photoionisation and photoheating          |
|                                     | |         processes just use the optically thin rates.                       |
|                                     | | ``1`` - Switches on self-shielding. The photoionisation and photoheating   |
|                                     | |         rates are attenuated according to the shielding length given in    |
|                                     | |         the **Gas Variables** structure (see below).                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``N_spectra``                       | | The number of UV spectra. You will need to provide a cross sections data   |
|                                     | | file for each spectrum (see the :ref:`ChimesData_label` section). Maximum  |
|                                     | | number of spectra: 20.                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``redshift_dependent_``           | | Gives the index of the UV spectrum for the redshift-dependent UV           |
| | ``UVB_index``                     | | background (starting from 0). If set to -1, no redshift-dependent UVB is   |
|                                     | | included.                                                                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``use_redshift_``                 | | If this integer flag is set to 1, the equilibrium tables are also          |
| | ``dependent_eqm_``                | | interpolated to the current redshift when using the redshift-dependent     |
| | ``tables``                        | | UVB, in addition to the cross sections tables (see the                     |
|                                     | | :ref:`RedshiftUVB_label` section).                                         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``redshift``                        | | The current redshift, used to interpolate the cross sections and           |
|                                     | | equilibrium tables when using the redshift-dependent UVB. If running a     |
|                                     | | cosmological simulation, this will need to be updated at each              |
|                                     | | hydrodynamics time-step, but is then held fixed while looping over all     |
|                                     | | active gas particles/cells.                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``reionisation_``                 | | Sets the redshift of reionisation. The redshift-dependent UVB, if being    |
| | ``redshift``                      | | used, is set to zero before reionisation.                                  |
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
|                                     | | and in the cooling rate due to energy transfer between gas and dust. In    |
|                                     | | CHIMES, we use a fixed grain temperature for all gas particles/cells;      |
|                                     | | there is currently no option to calculate this separately for each         |
|                                     | | particle/cell.                                                             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``cmb_temperature``                 | | Temperature of the Cosmic Microwave Background, used for Compton           |
|                                     | | cooling from the CMB. Note that CHIMES does not automatically set this     |
|                                     | | according to the specified redshift. If running a cosmological             |
|                                     | | simulation, it will need to be updated each hydrodynamical time-step.      |
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
| ``element_included[9]``             | | Array of integer flags (each ``0`` or ``1``) that control whether to       |
|                                     | | include each metal element in the chemical network. These are given in the |
|                                     | | order: C, N, O, Ne, Mg, Si, S, Ca, Fe. If a given element is not included, |
|                                     | | all of the ions and molecules involving that element, along with all       |
|                                     | | associated reactions, are removed from the network. Note that hydrogen     |
|                                     | | and helium are always included, so there are no flags for these.           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``speciesIndices[``               | | Array of indices that maps each species included in the network to its     |
| | ``CHIMES_TOTSIZE]``               | | position in the abundance array in the **Gas Variables** structure (see    |
|                                     | | below). This array is constructed by the ``init_chimes()`` routine when    |
|                                     | | you initialise CHIMES, based on the ``element_included`` flags that        |
|                                     | | determine which elements are included in the network, so they do not need  |
|                                     | | to be provided directly by the user. If a species is not included in the   |
|                                     | | network, its index is set to ``-1``. ``CHIMES_TOTSIZE`` is the total       |
|                                     | | number of species in the full network (i.e. ``157``).                      |
|                                     | | To find the position of a given species, you can use the enumerated        |
|                                     | | species names given at the end of ``chimes_proto.h``. This enumerates      |
|                                     | | where that species appears in the full network, and then it can be looked  |
|                                     | | up in the ``speciesIndices[]`` array to find its position in the reduced   |
|                                     | | network. For example, the following would give you the position of H2 in   |
|                                     | | the reduced network: ``speciesIndices[sp_H2]``. The names of all species   |
|                                     | | here have been prepended with ``sp_`` (for species) to reduce the risk     |
|                                     | | that a species name also appears as a variable name in the hydrodynamics   |
|                                     | | code.                                                                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``totalNumberOf``                 | | The total number of species in the reduced network. This is calculated by  |
| | ``Species``                       | | the ``init_chimes()`` routine when you initialise CHIMES, based on the     |
|                                     | | ``element_included`` flags that determine which elements are included in   |
|                                     | | the network, so it does not need to be provided directly by the user.      |
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
|                                     | | non-equilibrium network using pre-computed cooling tables. See the         |
|                                     | | :ref:`HybridCool_label` section for details on how to use this option.     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``*hybrid_data``                    | | A void pointer that can be used to point to the User-defined structure     |
|                                     | | containing additional data needed for the User's hybrid cooling function.  |
|                                     | | Only used if ``hybrid_cooling_mode == 1`` (see the :ref:`HybridCool_label` |
|                                     | | section for details).                                                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``*hybrid_cooling_fn``              | | A function pointer to the User-defined hybrid cooling function. Only used  |
|                                     | | if ``hybrid_cooling_mode == 1`` (see the :ref:`HybridCool_label` section   |
|                                     | | for details).                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``*allocate_gas_``                | | A function pointer to the User-defined function that allocates memory to   |
| | ``hybrid_data_fn``                | | the ``*hybrid_data`` structure in ``gasVariables`` (not the one in         |
|                                     | | ``globalVariables``). Only used if ``hybrid_cooling_mode == 1`` (see the   |
|                                     | | :ref:`HybridCool_label` section for details).                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``*free_gas_``                    | | A function pointer to the User-defined function that frees memory          |
| | ``hybrid_data_fn``                | | allocated for the ``*hybrid_data`` structure in ``gasVariables`` (not the  |
|                                     | | one in ``globalVariables``). Only used if ``hybrid_cooling_mode == 1``     |
|                                     | | (see the :ref:`HybridCool_label` section for details).                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

Gas Variables
^^^^^^^^^^^^^

The ``gasVariables`` structure contains variables specific to each gas particle/cell (e.g. density, temperature etc.). When you create an instance of this structure, you will need to use the ``allocate_gas_abundances_memory()`` routine from ``init_chimes.c`` to allocate memory for all of the dynamically allocated arrays in the structure, such as the abundance array. Then when you are finished with a given instance, you will need to use the ``free_gas_abundances_memory()`` routine to free that memory. 

There are two ways in which you could implement this structure in the hydrodynamics code. Firstly, you could create a separate instance of ``gasVariables`` for every gas particle/cell in the simulation. This would store each particle/cell's abundance array, and then every hydrodynamical time-step when you come to integrate the chemistry and cooling, you simply update the various gas variables in this structure from that particle/cell's hydrodynamic quantities and then call the CHIMES chemistry solver routine on that structure. 

The main downside of this approach is that the ``gasVariables`` structure contains dynamically allocated arrays, for example for the abundance array. These need to be allocated in this way because the size of the array is typically determined at run-time by the parameters provided by the User, for example to set which elements to include in the network. However, this makes it more complicated when moving a particle/cell from one MPI task to another, because you cannot just send the ``gasVariables`` struct for that particle/cell. Instead you have to read the dynamically allocated arrays into separate buffers, free the memory from those arrays, send the ``gasVariables`` struct and all of the array buffers to the new MPI task, allocate memory for those arrays on the new MPI task, and then read the data from the buffer back into the corresponding arrays in the ``gasVariables`` structure. 

Alternatively, you could just add an array of fixed size for the abundances to the existing hydrodynamics data structure for each gas particle/cell. Then every time-step when you call the CHIMES chemistry solver for each particle/cell, you first create an instance of the ``gasVariables`` structure, allocate memory for its arrays using the ``allocate_gas_abundances_memory()`` routine, copy over the abundances and all of the other hydrodynamic quantities and then call the chemistry solver on that structure. Then at the end of the time-step you can update that particle/cell's abundance array from the ``gasVariables`` structure, and then free the memory that was allocated to the arrays in ``gasVariables`` using the ``free_gas_abundances_memory()`` routine. 

This means that you would need to know the size of the abundance arrays at compile-time rather than run-time, so you would need to re-compile the code whenever you changed which elements to include in the network. It also means that you are freeing and allocating memory in the ``gasVariables`` structures a lot more, as you have to do it every time-step for every active gas particle/cell, although in practice the extra expense from doing this is very small compared to the cost of the chemistry solver itself. 

The variables in the ``gasVariables`` structure are described below. 

+-------------------------------------+------------------------------------------------------------------------------+
| Variable                            | Description                                                                  |
+=====================================+==============================================================================+
| | ``element_``                      | | Abundance of each element relative to hydrogen, i.e. ``n_i / n_Htot``,     |
| | ``abundances[10]``                | | where ``n_i`` and ``n_Htot`` are the number densities of element ``i`` and |
|                                     | | hydrogen, respectively. These are given in the order He, C, N, O, Ne, Mg,  |
|                                     | | Si, S, Ca, and Fe. Note that the abundance of hydrogen is unity by         |
|                                     | | definition, so it is not included in this array.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``nH_tot``                          | | Total hydrogen number density, in units of cm^-3.                          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``temperature``                     | | Gas temperature, in units of K.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TempFloor``                       | | Temperature floor used in the cooling integration, in K.                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``divVel``                          | | Divergence of the velocity field, in cgs units. Only used if               |
|                                     | | ``globalVariables.StaticMolCooling == 0``.                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``doppler_broad``                   | | Doppler broadening parameter due to turbulence, in km/s, used in the H2    |
|                                     | | self-shielding function (see Richings et al. 2014b).                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``*isotropic_photon_``            | | The strength of the ionising radiation field for each spectrum, defined as |
| | ``density``                       | | the number density of hydrogen-ionising photons (i.e. those with energies  |
|                                     | | >13.6 eV), in units of ``photons cm^-3``. The flux of ionising photons,    |
|                                     | | in units of ``photons cm^-2 s^-1``, can then be found by multiplying       |
|                                     | | ``isotropic_photon_density`` by the speed of light. This is a dynamically  |
|                                     | | allocated array of length ``globalVariables.N_spectra``.                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``*G0_parameter``                   | | The strength of the radiation field for each spectrum in the 6-13.6 eV     |
|                                     | | band, ``G0``, in Habing units, divided by the flux of ionising photons,    |
|                                     | | i.e.: ``G0_parameter = G0 / (isotropic_photon_density *``                  |
|                                     | | ``speed_of_light)``.                                                       |
|                                     | | By normalising the ``G0_parameter`` in this way, we can vary the           |
|                                     | | normalisation of the whole UV spectrum, in both the 6-13.6 eV and the      |
|                                     | | >13.6 eV bands, by just changed in the ``isotropic_photon_density``.       |
|                                     | | The ``G0_parameter`` then only depends on the shape of the spectrum.       |
|                                     | | This is a dynamically allocated array of length                            |
|                                     | | ``globalVariables.N_spectra``                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``*H2_dissocJ``                     | | The strength of the radiation field for each spectrum in the               |
|                                     | | 12.24-13.51 eV band, defined as the number density of photons in this band |
|                                     | | divided by the ``isotropic_photon_density`` parameter times the speed of   |
|                                     | | light. This is used to calculate the photodissociation rate of H2. This is |
|                                     | | a dynamically allocated array of length ``globalVariables.N_spectra``.     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``cr_rate``                         | | Ionisation rate of HI due to cosmic rays, in units of s^-1.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``metallicity``                     | | Gas-phase metallicity, given as ``Z / Zsol``, where ``Zsol = 0.0129`` is   |
|                                     | | the solar metallicity. This is used to interpolate the equilibrium         |
|                                     | | abundance tables when ``ForceEqOn == 1``.                                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``dust_ratio``                      | | Dust-to-gas ratio, ``D/G``, divided by the Milky Way Dust-to-gas ratio of  |
|                                     | | ``D/G_MW = 0.06``.                                                         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``cell_size``                       | | Shielding length of the gas particle/cell, in units of cm. If              |
|                                     | | ``globalVariables.cellSelfShieldingOn == 1``, the photoionisation and      |
|                                     | | photoheating rates are attenuated by column densities of HI, H2, HeI, HeII |
|                                     | | and CO given by multiplying the number densities of each species by the    |
|                                     | | shielding length (see Richings et al. 2014b).                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``hydro_timestep``                  | | Total time over which to integrate the chemistry and cooling for each call |
|                                     | | to ``chimes_network()``, in seconds. This would be set to the              |
|                                     | | hydrodynamical time-step for the given gas particle/cell from the          |
|                                     | | hydrodynamics code. The CVODE solver will then sub-cycle this time-step to |
|                                     | | integrate the chemistry and cooling.                                       |
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
|                                     | | ``0`` - If the temperature is less than ``TempFloor``, the rate of change  |
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
| ``InitIonState``                    | | Integer that controls the initial chemical state that the abundance array  |
|                                     | | is set to when we call the ``initialise_gas_abundances()`` routine from    |
|                                     | | ``init_chimes.c``. This routine is a convenient way to set the abundances  |
|                                     | | to some initial state before we start the chemistry integration. Each      |
|                                     | | element is then set to the ionisation state given by this parameter. For   |
|                                     | | example, if set to ``0`` all elements will be neutral, whereas if set to   |
|                                     | | ``1`` all elements will be singly ionised.                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``constant_heating_``             | | Constant heating rate (positive for heating, negative for cooling), in     |
| | ``rate``                          | | units of erg cm^-3 s^-1, that is added on to the heating and cooling       |
|                                     | | rates from the CHIMES network. This can be used, for example, to include   |
|                                     | | the ``p dV`` work from adiabatic expansion/contraction in the cooling      |
|                                     | | integration. Set this to zero to just use the heating and cooling rates    |
|                                     | | from the network.                                                          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``*abundances``                     | | Array containing the abundances of each species in the network, given as   |
|                                     | | ``n_i / n_Htot`` where ``n_i`` and ``n_Htot`` are the number densities of  |
|                                     | | species ``i`` and hydrogen, respectively. This is a dynamically allocated  |
|                                     | | array of length ``globalVariables.totalNumberOfSpecies``.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``*hybrid_data``                    | | A void pointer that can be used to point to the User-defined structure     |
|                                     | | containing additional data needed for the User's hybrid cooling function.  |
|                                     | | This is separate from the ``*hybrid_data`` in the **globalVariables**      |
|                                     | | structure, and is used for additional variables that are specific to each  |
|                                     | | gas particle/cell. This is only used if ``hybrid_cooling_mode == 1`` (see  |
|                                     | | the :ref:`HybridCool_label` section for details).                          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

References
^^^^^^^^^^

| `Richings et al. (2014b) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.2780R>`_

