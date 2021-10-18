.. CHIMES cross sections 
   Alexander Richings, 3rd March 2020 

.. _ChimesCrossSections_label: 

cross_sections
--------------

For each UV spectrum that you want to include in the chemical modelling, you will need to specify a cross sections data file. This contains the various photoionisation cross sections, shielding factors, and photoheating rates for each photo-chemical reaction, averaged over the given spectrum. Quantities corresponding to individual reactions are organised into the same reaction groups as used in ``chimes_main_data.hdf5`` (see the :ref:`ChimesMainData_label` section for more details). Note that CHIMES assumes that the reactions in each group are arranged in the same order in ``chimes_main_data.hdf5`` and the individual cross sections data files. 

The ``chimes-data`` repository includes several cross sections data files for commonly used spectra, including the Black (1987) local interstellar radiation field (``cross_Sections_B87.hdf5``), the Sazonov et al (2004) average quasar spectrum (``cross_sections_S04.hdf5``), the redshift-dependent extragalactic backgrounds from Haardt & Madau (2001) (``HM01_cross_sections/``), Haardt & Madau (2012) (``HM12_cross_sections/``) and Faucher-Giguere (2020) (``FG19_cross_sections/``). We also include several UV spectra from a single stellar population using Starburst 99 models (Leitherer et al. 2014) at different stellar ages (``starburstCrossSections/``; see Richings et al. in prep for details). 

If you want to include your own UV spectra in CHIMES, you will need to create the corresponding cross sections file. We are currently putting together a python script to create this file for any user-defined spectrum; we will make this availabe in due course. 

The data arrays included in the cross sections files are described in more detail below. 

photoion_fuv
^^^^^^^^^^^^

Photoionisation of species with an ionisation energy <13.6 eV. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. The first number gives the number of reactions that do not     |
|                                     | | involve molecules, while the second number gives the total number of       |
|                                     | | reactions, including molecules.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``sigmaPhot``                       | | A 1-dimensional array of length ``N_reactions`` that gives the             |
|                                     | | photoionisation cross section of each reaction, averaged over the given    |
|                                     | | UV spectrum (see equation 2.5 in Richings et al. 2014a).                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``epsilonPhot``                     | | A 1-dimensional array of length ``N_reactions`` that gives the spectrum-   |
|                                     | | averaged excess energy of ionising photons in the optically thin limit,    |
|                                     | | used to calculate the photoheating rate (see equation 3.6 in Richings et   |
|                                     | | al. 2014a).                                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

photoion_euv
^^^^^^^^^^^^

Photoionisation of species with an ionisation energy >13.6 eV. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. The first number gives the number of reactions that do not     |
|                                     | | involve molecules, while the second number gives the total number of       |
|                                     | | reactions, including molecules.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``sigmaPhot``                       | | A 1-dimensional array of length ``N_reactions`` that gives the             |
|                                     | | photoionisation cross section of each reaction, averaged over the given    |
|                                     | | UV spectrum (see equation 2.5 in Richings et al. 2014a).                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``shieldFactor_1D``                 | | A 3-dimensional array of length (``N_reactions`` x 3 x                     |
|                                     | | ``N_Column_densities``). In the second dimension of this array, index 0    |
|                                     | | contains the factor ``S_gas1`` from equation 3.8 of Richings et al.        |
|                                     | | (2014b) for each reaction as a function of HI column density, for          |
|                                     | | calculating the shielded photoionisation rates; while indices 1 and 2      |
|                                     | | contain the first integrals in the numerator and denominator,              |
|                                     | | respectively, of equation 3.22 in Richings et al. (2014b) for each         |
|                                     | | reaction as a function of HI column density, for calculating the shielded  |
|                                     | | photoheating rates.                                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``shieldFactor_2D``                 | | A 4-dimensional array of length (``N_reactions`` x 6 x                     |
|                                     | | ``N_Column_densities`` x ``N_Column_densities``). In the second dimension  |
|                                     | | of this array, indices 0 and 1 contain the factors ``S_gas2`` and          |
|                                     | | ``S_gas3`` from equation 3.8 of Richings et al. (2014b) for each reaction, |
|                                     | | for calculating the shielded photoionisation rates. Indices 2-5 contain    |
|                                     | | the second and third integrals in the numerator and denominator of         |
|                                     | | equation 3.22 in Richings et al. (2014b) for each reaction, for            |
|                                     | | calculating the shielded photoheating rates.                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

photoion_auger_fuv
^^^^^^^^^^^^^^^^^^

Auger photoionisation of species with an ionisation energy <13.6 eV. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. The first number gives the number of reactions that do not     |
|                                     | | involve molecules, while the second number gives the total number of       |
|                                     | | reactions, including molecules.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``sigmaPhot``                       | | A 1-dimensional array of length ``N_reactions`` that gives the             |
|                                     | | photoionisation cross section of each reaction, averaged over the given    |
|                                     | | UV spectrum (see equation 2.5 in Richings et al. 2014a).                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

photoion_auger_euv
^^^^^^^^^^^^^^^^^^

Auger photoionisation of species with an ionisation energy >13.6 eV. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. The first number gives the number of reactions that do not     |
|                                     | | involve molecules, while the second number gives the total number of       |
|                                     | | reactions, including molecules.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``sigmaPhot``                       | | A 1-dimensional array of length ``N_reactions`` that gives the             |
|                                     | | photoionisation cross section of each reaction, averaged over the given    |
|                                     | | UV spectrum (see equation 2.5 in Richings et al. 2014a).                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

General Spectrum Data
^^^^^^^^^^^^^^^^^^^^^

The following data arrays contain general information about the given UV spectrum, and are found in the root of the cross sections HDF5 data file. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``isotropic_photon_density``        | | The strength of the ionising radiation field, defined as the number        |
|                                     | | density of hydrogen-ionising photons (i.e. those with energies >13.6 eV),  |
|                                     | | in units of ``cm^-3``. The flux of ionising photons, in units of           |
|                                     | | ``photons cm^-2 s^-1``, can then be found by multiplying                   |
|                                     | | ``isotropic_photon_density`` by the speed of light. We have included       |
|                                     | | this  parameter in the data files for convenience. However, the user is    |
|                                     | | free to set the corresponding ``isotropic_photon_density`` parameter in    |
|                                     | | the **gasVariables** structure to a different value if they wish, for      |
|                                     | | example to vary the strength of the radiation field compared to the        |
|                                     | | fiducial value given in the data files.                                    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``G0_parameter``                    | | The strength of the radiation field in the 6-13.6 eV band, ``G0``, in      |
|                                     | | Habing units, divided by the flux of ionising photons, i.e.:               |
|                                     | | ``G0_parameter = G0 / (isotropic_photon_density *``                        |
|                                     | | ``speed_of_light)``.                                                       |
|                                     | | By normalising the ``G0_parameter`` in this way, we can vary the           |
|                                     | | normalisation of the whole UV spectrum, in both the 6-13.6 eV and the      |
|                                     | | >13.6 eV bands, by just changed in the ``isotropic_photon_density``.       |
|                                     | | The ``G0_parameter`` then only depends on the shape of the spectrum.       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``H2_dissocJ``                      | | The strength of the radiation field in the 12.24-13.51 eV band, defined as |
|                                     | | the number density of photons in this band divided by the                  |
|                                     | | ``isotropic_photon_density`` parameter times the speed of light. This is   |
|                                     | | used to calculate the photodissociation rate of H2.                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

TableBins
^^^^^^^^^

The **TableBins** group contains the Column Density table bins that are used to tabulate and interpolate the various shielding factor arrays. 

References
^^^^^^^^^^

| `Black (1987) <https://ui.adsabs.harvard.edu/abs/1987ASSL..134..731B>`_ 
| `Faucher-Giguere (2020) <https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1614F>`_ 
| `Haardt & Madau (2001) <https://ui.adsabs.harvard.edu/abs/2001cghr.confE..64H>`_ 
| `Haardt & Madau (2012) <https://ui.adsabs.harvard.edu/abs/2012ApJ...746..125H>`_ 
| `Leitherer et al. (2014) <https://ui.adsabs.harvard.edu/abs/2014ApJS..212...14L>`_ 
| `Richings et al. (2014a) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.440.3349R>`_
| `Richings et al. (2014b) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.2780R>`_
| `Sazonov et al. (2004) <https://ui.adsabs.harvard.edu/abs/2004MNRAS.347..144S>`_ 
