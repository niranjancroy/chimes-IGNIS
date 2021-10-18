.. CHIMES Driver UV Field
   Alexander Richings, 5th March 2020 

.. _ChimesDriverUVField_label:

UV Field
--------

The UV radiation field used in the chemistry solver can be specified via the ``UV_field`` parameter. The possible options are described in detail below. 

+-------------------------------------+------------------------------------------------------------------------------+
| UV Field                            | Description                                                                  |
+=====================================+==============================================================================+
| ``None``                            | | No UV radiation field is used.                                             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``HM01``                            | | Uses the extragalactic UV background from Haardt & Madau (2001). Note      |
|                                     | | that ``chimes-data`` currently only includes the cross-sections data at    |
|                                     | | redshift zero for this UVB. The normalisation of this radiation field      |
|                                     | | from the cross sections table is multiplied by the                         |
|                                     | | ``radiation_field_normalisation_factor`` from the parameter file.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``HM12``                            | | Uses the extragalactic UV background from Haardt & Madau (2012). It will   |
|                                     | | interpolate the redshift-dependent cross sections and radation field       |
|                                     | | strength from the HM12 cross sections tables in ``chimes-data`` to the     |
|                                     | | current redshift as given by the ``redshift`` parameter. If                |
|                                     | | ``redshift > reionisation_redshift``, the UVB is set to zero.              |
|                                     | | The normalisation of this radiation field from the cross sections table    |
|                                     | | is multiplied by the ``radiation_field_normalisation_factor`` from         |
|                                     | | the parameter file.                                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``FG20``                            | | Uses the extragalactic UV background from Faucher-Giguere (2020). It will  |
|                                     | | interpolate the redshift-dependent cross sections and radation field       |
|                                     | | strength from the FG20 cross sections tables in ``chimes-data`` to the     |
|                                     | | current redshift as given by the ``redshift`` parameter. If                |
|                                     | | ``redshift > reionisation_redshift``, the UVB is set to zero.              |
|                                     | | The normalisation of this radiation field from the cross sections table    |
|                                     | | is multiplied by the ``radiation_field_normalisation_factor`` from         |
|                                     | | the parameter file.                                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``B87``                             | | Uses the interstellar radiation field in the local solar neighbourhood in  |
|                                     | | the Milky Way, from Black (1987). The normalisation of this radiation      |
|                                     | | field from the cross sections table is multiplied by the                   |
|                                     | | ``radiation_field_normalisation_factor`` from the parameter file.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``Colibre``                         | | Uses the SP20 redshift-dependent extragalactic UV background plus the      |
|                                     | | interstellar radiation field (ISRF) from Black (1987) scaled relative to   |
|                                     | | ISRF in the Milky Way based on the Jeans column density (see Ploeckinger   |
|                                     | | & Schaye 2020). Note that the ``radiation_field_normalisation_factor``     |
|                                     | | parameter is not used here. Instead, the ISRF component can be             |
|                                     | | re-normalised using the ``colibre_scale_MW_ISRF`` parameter. Ploeckinger   |
|                                     | | & Schaye (2020) use a fiducial value of ``0.1`` for this parameter.        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``StellarFluxes``                   | | Uses the FG20 extragalactic UV background at redshift zero, plus the       |
|                                     | | stellar fluxes from the star particles using the UV spectra from           |
|                                     | | Starburst 99 models (Leitherer et al. 2014) in 8 stellar age bins. See     |
|                                     | | Richings et al. (in prep) for details. Note that the                       |
|                                     | | ``radiation_field_normalisation_factor`` parameter is not used here.       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``S04``                             | | Uses the average quasar UV spectrum from Sazonov et al. (2004). The        |
|                                     | | normalisation of the spectrum is determined by the bolometric AGN          |
|                                     | | luminosity, specified via the ``bolometric_AGN_luminosity_cgs``            |
|                                     | | parameter, and the distance to the AGN. If ``IO_mode`` is set to           |
|                                     | | ``grid``, the distance is specified via the ``distance_to_AGN_kpc``        |
|                                     | | parameter. If ``IO_mode`` is set to ``snapshot``, we instead specify       |
|                                     | | the position of the AGN via the ``AGN_position_x_kpc``,                    |
|                                     | | ``AGN_position_y_kpc`` and ``AGN_position_z_kpc`` parameters. The          |
|                                     | | distance from the AGN to each gas particle is then calculated separately.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

References
^^^^^^^^^^

| `Black (1987) <https://ui.adsabs.harvard.edu/abs/1987ASSL..134..731B>`_ 
| `Faucher-Giguere (2020) <https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1614F>`_ 
| `Haardt & Madau (2001) <https://ui.adsabs.harvard.edu/abs/2001cghr.confE..64H>`_ 
| `Haardt & Madau (2012) <https://ui.adsabs.harvard.edu/abs/2012ApJ...746..125H>`_ 
| `Leitherer et al. (2014) <https://ui.adsabs.harvard.edu/abs/2014ApJS..212...14L>`_ 
| `Ploeckinger & Schaye (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200614322P/abstract>`_
| `Sazonov et al. (2004) <https://ui.adsabs.harvard.edu/abs/2004MNRAS.347..144S/abstract>`_

