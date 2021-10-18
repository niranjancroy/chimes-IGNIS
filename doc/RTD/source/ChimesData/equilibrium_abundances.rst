.. CHIMES equilibrium abundances 
   Alexander Richings, 3rd March 2020 

.. _ChimesEqmAbundances_label: 

equilibrium_abundances
----------------------

The equilibrium abundances data files tabulate the equilibrium chemical abundances as a function of temperature, density and metallicity. These are used when the ``gasVariables.ForceEqOn`` parameter is set to ``1``. 

We have included several examples of these pre-computed equilibrium tables, for different model assumptions, including for the COLIBRE model for the radiation field, shielding and dust depletion (Ploeckinger & Schaye 2020), both with only hydrogen and helium in the network (``EqAbundancesTables/colibre_HHe``) and with the full network but only at high metallicities (``EqAbundancesTables/colibre_iso_full.hdf5``), for the redshift zero extragalactic UV background from Haardt & Madau (2001), but with and without shielding (``EqAbundancesTables/HM01/EqAbundances_HM01.hdf5`` and ``EqAbundancesTables/HM01/EqAbundances_HM01_noShield.hdf5``, respectively), and for the redshift zero extragalactic UV background from Haardt & Madau (2012) without shielding (``EqAbundancesTables/HM12/EqAbundances_HM12_noShield.hdf5``).

Note that, if you do load an equilibrium abundance table, the number of species in the table needs to match the number of species that are included in the CHIMES network, for example if you have switched off individual elements. 

New equilibrium abundance tables can be created using the ``chimes-driver.py`` script. See the :ref:`ChimesDriver_label` section for more details. 

The data arrays included in the equilibrium abundance tables are described in detail below. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``Abundances``                      | | A 4-dimensional array of size (``N_Temperatures`` x ``N_Densities`` x      |
|                                     | | ``N_Metallicities`` x ``N_species``), giving the abundance of each         |
|                                     | | species as a function of temperature, density and metallicity. The         |
|                                     | | abundances are stored as ``log10(n_i / n_elem)``, where ``n_i`` is the     |
|                                     | | number density of the ion or molecule, and ``n_elem`` is the total         |
|                                     | | number density of the corresponding element.                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TableBins/N_Temperatures``        | | The number of temperature bins.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TableBins/N_Densities``           | | The number of density bins.                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TableBins/N_Metallicities``       | | The number of metallicity bins.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TableBins/N_species``             | | The number of species included in the network.                             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TableBins/Temperatures``          | | A 1-dimensional array of length ``N_Temperatures`` containing the          |
|                                     | | temperature bins, given as log10(T [K]).                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TableBins/Densities``             | | A 1-dimensional array of length ``N_Densities`` containing the density     |
|                                     | | bins, given as log10(nHtot [cm^-3]).                                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``TableBins/Metallicities``         | | A 1-dimensional array of length ``N_Metallicities`` containing the         |
|                                     | | metallicity bins, given as log10(Z / Zsol), where Zsol = 0.0129 is the     |
|                                     | | solar metallicity.                                                         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

References
^^^^^^^^^^
 
| `Haardt & Madau (2001) <https://ui.adsabs.harvard.edu/abs/2001cghr.confE..64H>`_ 
| `Haardt & Madau (2012) <https://ui.adsabs.harvard.edu/abs/2012ApJ...746..125H>`_ 
| `Ploeckinger & Schaye (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200614322P/abstract>`_
