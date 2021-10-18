.. CHIMES Dust Depletion
   Alexander Richings, 6th March 2020

.. _ChimesDriverDustDepletion_label:

Dust Depletion
--------------

The ``dust_depletion`` parameter controls how we model the depletion of metals on to dust grains in the chemistry solver. The possible options are described below. 

+-------------------------------------+------------------------------------------------------------------------------+
| Dust Depletion Model                | Description                                                                  |
+=====================================+==============================================================================+
| ``None``                            | | The dust-to-gas ratio, is scaled linearly with the metallicity, i.e. we    |
|                                     | | assume a constant dust-to-metals ratio. However, the abundances of metals  |
|                                     | | in the gas phase are not reduced to account for the depletion of metals    |
|                                     | | on to dust grains.                                                         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``J09``                             | | The depletion of individual metals on to dust grains is calculated from    |
|                                     | | the observationally measured depletion factors (i.e. the fraction of each  |
|                                     | | element that is in dust) in the Milky Way from Jenkins (2009). These are   |
|                                     | | expressed as functions of a parameter F_star, which represents the overal  |
|                                     | | strength of dust depletion and is itself a function of density,            |
|                                     | | ``n_Htot``. We also impose an upper temperature cut of 1e6 K, above which  |
|                                     | | all metals are put in the gas phase, to account for the destruction of     |
|                                     | | dust grains due to sputtering at high gas temperatures. The total dust-    |
|                                     | | to-gas ratio is then calculated by summing over the depletions of each     |
|                                     | | individual element. The dust abundances and depletion factors are thus     |
|                                     | | functions of gas temperature, density and metallicity. See Richings et     |
|                                     | | al. (in prep) for details of how we implement this model in CHIMES.        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``DC16``                            | | As ``J09``, but using the updated fit parameters from De Cia et al. (2016) |
|                                     | | for the individual element depletion factors, where available. For         |
|                                     | | elements that are not included in De Cia et al. (2016), we use the fits    |
|                                     | | from Jenkins (2009). See Richings et al. (in prep) for details.            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``Colibre``                         | | The dust-to-gas ratio and depletion factors of individual elements are     |
|                                     | | calculated as in the COLIBRE model. See Ploeckinger & Schaye (2020) for    |
|                                     | | details.                                                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

References
^^^^^^^^^^

| `De Cia et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016A%26A...596A..97D>`_
| `Jenkins (2009) <https://ui.adsabs.harvard.edu/abs/2009ApJ...700.1299J>`_
| `Ploeckinger & Schaye (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200614322P/abstract>`_
