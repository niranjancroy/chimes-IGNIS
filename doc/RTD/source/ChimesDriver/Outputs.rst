.. CHIMES Driver Outputs
   Alexander Richings, 4th March 2020 

.. _ChimesDriverOutputs_label:

Outputs
-------

The outputs that are produced by CHIMES Driver depend on the combination of the ``IO_mode`` and ``driver_mode`` parameters. The possible variations are described below. Note that, if the output file, as specified by the ``output_file`` parameter, is already present, then CHIMES Driver will check at the start whether the output arrays are already present in this file. If they are, it will exit with an error message, rather than try to over-write these arrays. 

Snapshot
^^^^^^^^ 

If ``IO_mode == snapshot``, the ``output_file`` can be set to the same as the ``input_file``, in which case the output arrays will be appended to the input snapshot, or it can be set to a new file. The arrays that are written out to the output file depend on the ``driver_mode`` as follows: 

+-------------------------------------+------------------------------------------------------------------------------+
| ``driver_mode``                     | Output Arrays                                                                |
+=====================================+==============================================================================+
| ``eqm_state``                       | | ``<hdf5_output_group>/EqmChemistryAbundances`` - A 2-dimensional array     |
|                                     | |    of size (``N_gas`` x ``N_species``), where ``N_gas`` is the number of   |
|                                     | |    gas particles and ``N_species`` is the number of species in the         |
|                                     | |    network, that gives the final equilibrium abundances, relative to       |
|                                     | |    hydrogen, for each gas particle.                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``cooling_rates``                   | | ``<hdf5_output_group>/log_cooling_rate`` - A 1-dimensional array of        |
|                                     | |    length ``N_gas``, where ``N_gas`` is the number of gas particles, that  |
|                                     | |    gives the log10 of the cooling rate (i.e. summed over all cooling       |
|                                     | |    channels) of each particle for the final abundances, in units of        |
|                                     | |    erg cm^-3 s^-1.                                                         |
|                                     | | ``<hdf5_output_group>/log_heating_rate`` - A 1-dimensional array of        |
|                                     | |    length ``N_gas``, where ``N_gas`` is the number of gas particles, that  |
|                                     | |    gives the log10 of the heating rate (i.e. summed over all heating       |
|                                     | |    channels) of each particle for the final abundances, in units of        |
|                                     | |    erg cm^-3 s^-1.                                                         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``noneq_evolution``                 | | ``<hdf5_output_group>/AbundanceEvolution`` - A 3-dimensional array         |
|                                     | |    of size (``N_gas`` x ``N_species`` x ``N_time``), where ``N_gas`` is    |
|                                     | |    the number of gas particles, ``N_species`` is the number of species in  |
|                                     | |    the network, and ``N_time`` is the number of time outputs. This gives   |
|                                     | |    the abundances relative to hydrogen for each gas particle at each time  |
|                                     | |    output.                                                                 |
|                                     | | ``<hdf5_output_group>/TemperatureEvolution`` - A 2-dimensional array       |
|                                     | |    of size (``N_gas`` x ``N_time``), that gives the temperature of each    |
|                                     | |    gas particle at each output time.                                       |
|                                     | | ``<hdf5_output_group>/TimeArray_seconds`` - A 1-dimensional array of       |
|                                     | |    length ``N_time`` that gives the time of each output in seconds.        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

If ``UV_field == StellarFluxes`` and ``compute_stellar_fluxes == 1``, we also write out the stellar fluxes from the ``ChimesFluxIon_arr`` and ``ChimesFluxG0_arr`` arrays (see the :ref:`ChimesDriverInputs_label` section for how these arrays are defined). The names of the datasets that these arrays are written out to in the output HDF5 file are given by the ``snapshot_flux_ion_array`` and ``snapshot_flux_G0_array`` parameters (see the :ref:`ChimesDriverParam_label` section for details). If these parameters are set to ``None``, then the corresponding arrays are not written out. 

Grid
^^^^

If ``IO_mode == grid``, the output arrays will be written to ``output_file`` depending on the ``driver_mode`` as follows: 

+-------------------------------------+------------------------------------------------------------------------------+
| ``driver_mode``                     | Output Arrays                                                                |
+=====================================+==============================================================================+
| ``eqm_state``                       | | ``Abundances`` - A 4-dimensional array of size (``N_Temperatures`` x       |
|                                     | |    ``N_Densities`` x ``N_Metallicities`` x ``N_species``), where           |
|                                     | |    ``N_Temperatures``, ``N_Densities`` and ``N_Metallicities`` are the     |
|                                     | |    number of temperature, density and metallicity bins in the grid, and    |
|                                     | |    ``N_species`` is the number of species in the network. This gives the   |
|                                     | |    final equilibrium abundances relative hydrogen, stored as linear not    |
|                                     | |    log.                                                                    |
|                                     | | ``TableBins/N_Temperatures`` - Number of temperature bins in the grid.     |
|                                     | | ``TableBins/N_Densities`` - Number of density bins in the grid.            |
|                                     | | ``TableBins/N_Metallicities`` - Number of metallicity bins in the grid.    |
|                                     | | ``TableBins/Temperatures`` - A 1-dimensional array of length               |
|                                     | |    ``N_Temperatures`` containing the temperature bins, given as            |
|                                     | |    log10(T [K]).                                                           |
|                                     | | ``TableBins/Densities`` - A 1-dimensional array of length ``N_Densities``  |
|                                     | |    containing the density bins, given as log10(nHtot [cm^-3]).             |
|                                     | | ``TableBins/Metallicities`` - A 1-dimensional array of length              |
|                                     | |    ``N_Metallicities`` containing the metallicity bins, given as           |
|                                     | |    log10(Z / Zsol), where Zsol = 0.0129 is the solar metallicity.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``eqm_table``                       | | ``Abundances`` - A 4-dimensional array of size (``N_Temperatures`` x       |
|                                     | |    ``N_Densities`` x ``N_Metallicities`` x ``N_species``), where           |
|                                     | |    ``N_Temperatures``, ``N_Densities`` and ``N_Metallicities`` are the     |
|                                     | |    number of temperature, density and metallicity bins in the grid, and    |
|                                     | |    ``N_species`` is the number of species in the network. This gives the   |
|                                     | |    final equilibrium abundances relative to the abundance of the           |
|                                     | |    corresponding element, stored as log10(n_i / n_elem). By recording the  |
|                                     | |    abundances in this way, the resulting table can be used as an input     |
|                                     | |    equilibrium table to the CHIMES solver, via the                         |
|                                     | |    ``EqAbundanceTable_filename`` parameter.                                |
|                                     | | ``TableBins/N_Temperatures`` - Number of temperature bins in the grid.     |
|                                     | | ``TableBins/N_Densities`` - Number of density bins in the grid.            |
|                                     | | ``TableBins/N_Metallicities`` - Number of metallicity bins in the grid.    |
|                                     | | ``TableBins/Temperatures`` - A 1-dimensional array of length               |
|                                     | |    ``N_Temperatures`` containing the temperature bins, given as            |
|                                     | |    log10(T [K]).                                                           |
|                                     | | ``TableBins/Densities`` - A 1-dimensional array of length ``N_Densities``  |
|                                     | |    containing the density bins, given as log10(nHtot [cm^-3]).             |
|                                     | | ``TableBins/Metallicities`` - A 1-dimensional array of length              |
|                                     | |    ``N_Metallicities`` containing the metallicity bins, given as           |
|                                     | |    log10(Z / Zsol), where Zsol = 0.0129 is the solar metallicity.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``cooling_rates``                   | | ``log_cooling_rate`` - A 3-dimensional array of size (``N_Temperatures``   |
|                                     | |    x ``N_Densities`` x ``N_Metallicities``), where ``N_Temperatures``,     |
|                                     | |    ``N_Densities`` and ``N_Metallicities`` are the number of temperature,  |
|                                     | |    density and metallicity bins in the grid. This gives the total cooling  |
|                                     | |    rate (i.e. summed over all cooling channels) as                         |
|                                     | |    log10(rate [erg cm^-3 s^-1]).                                           |
|                                     | | ``log_heating_rate`` - A 3-dimensional array of size (``N_Temperatures``   |
|                                     | |    x ``N_Densities`` x ``N_Metallicities``), where ``N_Temperatures``,     |
|                                     | |    ``N_Densities`` and ``N_Metallicities`` are the number of temperature,  |
|                                     | |    density and metallicity bins in the grid. This gives the total heating  |
|                                     | |    rate (i.e. summed over all heating channels) as                         |
|                                     | |    log10(rate [erg cm^-3 s^-1]).                                           |
|                                     | | ``TableBins/N_Temperatures`` - Number of temperature bins in the grid.     |
|                                     | | ``TableBins/N_Densities`` - Number of density bins in the grid.            |
|                                     | | ``TableBins/N_Metallicities`` - Number of metallicity bins in the grid.    |
|                                     | | ``TableBins/Temperatures`` - A 1-dimensional array of length               |
|                                     | |    ``N_Temperatures`` containing the temperature bins, given as            |
|                                     | |    log10(T [K]).                                                           |
|                                     | | ``TableBins/Densities`` - A 1-dimensional array of length ``N_Densities``  |
|                                     | |    containing the density bins, given as log10(nHtot [cm^-3]).             |
|                                     | | ``TableBins/Metallicities`` - A 1-dimensional array of length              |
|                                     | |    ``N_Metallicities`` containing the metallicity bins, given as           |
|                                     | |    log10(Z / Zsol), where Zsol = 0.0129 is the solar metallicity.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``noneq_evolution``                 | | ``AbundanceEvolution`` - A 5-dimensional array of size                     |
|                                     | |    (``N_Temperatures`` x ``N_Densities`` x ``N_Metallicities``             |
|                                     | |    x ``N_species`` x ``N_time``), where ``N_Temperatures``,                |
|                                     | |    ``N_Densities`` and ``N_Metallicities`` are the number of temperature,  |
|                                     | |    density and metallicity bins in the grid, ``N_species`` is the number   |
|                                     | |    of species in the network, and ``N_time`` is the number of time         |
|                                     | |    outputs. This gives the abundances relative to hydrogen at each time    |
|                                     | |    output.                                                                 |
|                                     | | ``TemperatureEvolution`` - A 4-dimensional array of size                   |
|                                     | |    (``N_Temperatures`` x ``N_Densities`` x ``N_Metallicities``             |
|                                     | |    x ``N_time``), that gives the temperatures at each time output.         |
|                                     | | ``TimeArray_seconds`` - A 1-dimensional array of length ``N_time`` that    |
|                                     | |    gives the time of each output in seconds.                               |
|                                     | | ``TableBins/N_Temperatures`` - Number of temperature bins in the grid.     |
|                                     | | ``TableBins/N_Densities`` - Number of density bins in the grid.            |
|                                     | | ``TableBins/N_Metallicities`` - Number of metallicity bins in the grid.    |
|                                     | | ``TableBins/Temperatures`` - A 1-dimensional array of length               |
|                                     | |    ``N_Temperatures`` containing the temperature bins, given as            |
|                                     | |    log10(T [K]).                                                           |
|                                     | | ``TableBins/Densities`` - A 1-dimensional array of length ``N_Densities``  |
|                                     | |    containing the density bins, given as log10(nHtot [cm^-3]).             |
|                                     | | ``TableBins/Metallicities`` - A 1-dimensional array of length              |
|                                     | |    ``N_Metallicities`` containing the metallicity bins, given as           |
|                                     | |    log10(Z / Zsol), where Zsol = 0.0129 is the solar metallicity.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
