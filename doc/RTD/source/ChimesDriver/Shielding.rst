.. CHIMES Driver Shielding
   Alexander Richings, 6th March 2020

.. _ChimesDriverShielding_label:

Shielding
---------

The radiation field defined through the ``UV_field`` parameter represents the radiation that is incident on the outside of the gas cloud. However, we also need to consider the self-shielding of the gas cloud, which attenuates the incident radiation field. We do this by specifying a shielding length, ``L_sh``, which characterises the size of the cloud. This then defines the total hydrogen column density that attenuates the radiation field as ``N_Htot = n_Htot * L_sh``, where ``n_Htot`` is the total number density of hydrogen nuclei, and the column densities of individual species, e.g. ``N_HI = n_HI * L_sh`` etc., which depend on the current abundances of those species. The cross sections data files used in CHIMES then tabulate how the various photoionisation and photoheating rates are attenuated as a function of these column densities, which accounts for how the shape of the UV spectra changes with column density (i.e. spectral hardening) as well as the overall reduction in the strength of the radiation field. See Section 3 of Richings et al. (2014b) for details. 

In CHIMES Driver, we use the ``shield_mode`` parameter to specify how to calculate the shielding length that is then passed to the chemistry solver. The possible options are described below. 

+-------------------------------------+------------------------------------------------------------------------------+
| ``shield_mode``                     | Description                                                                  |
+=====================================+==============================================================================+
| ``None``                            | | Shielding is disabled, and the optically thin rates are used.              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``read-in``                         | | The shielding column density of each gas particle is read in from the      |
|                                     | | snapshot. This needs to be given as the total hydrogen column density,     |
|                                     | | ``N_Htot``, in units of cm^-2. The shielding length is then calculated as  |
|                                     | | ``L_sh = N_Htot / n_Htot``, where ``n_Htot`` is the total number density   |
|                                     | | of hydrogen nuclei. This option can only be used if                        |
|                                     | | ``IO_mode == snapshot``.                                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``Jeans``                           | | The shielding length of each particle, or each point in the grid, is set   |
|                                     | | to the Jeans length, as a function of gas temperature and density. This    |
|                                     | | is multiplied by the ``shield_length_factor`` parameter, and is limited    |
|                                     | | according to the ``max_shield_length`` parameter (in cm).                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``Colibre``                         | | The shielding length of each particle, or each point in the grid, is set   |
|                                     | | according to the COLIBRE model (see Ploeckinger & Schaye 2020 for          |
|                                     | | details).                                                                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

References
^^^^^^^^^^

| `Ploeckinger & Schaye (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200614322P/abstract>`_ 
| `Richings et al. (2014b) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.2780R>`_
