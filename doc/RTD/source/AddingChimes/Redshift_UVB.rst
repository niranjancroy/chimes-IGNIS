.. Redshift UVB
   Alexander Richings, 23rd March 2020

.. _RedshiftUVB_label:

Redshift-Dependent UVB
----------------------

We have included routines in CHIMES that can be used to interpolate a redshift-dependent extragalactic UV background (e.g. Haardt & Madau 2012, Faucher-Giguere 2020) to the current redshift. We describe below how to use these routines when running CHIMES in your own hydrodynamics code.

Global Variables
^^^^^^^^^^^^^^^^

Some of the variables in the ``globalVariables`` structure relate to the redshift-dependent UVB option. We discuss these below (see also the :ref:`DataStructures_label` section).

* ``redshift_dependent_UVB_index`` - In CHIMES, you can define multiple UV spectra. This index gives the position of the spectrum corresponding to the redshift-dependent UVB, for example in the arrays that store the cross sections for each spectrum. When you specify the paths to the cross sections data files, they need to be specified in the same order. For example, if you set this index to ``0``, then the first path that you give will need to be for the redshift-dependent UVB. If this index is set to be negative, then the redshit-dependent UVB is not used.

* ``PhotoIonTablePath[CHIMES_MAX_UV_SPECTRA][500]`` - When you give the path to the cross sections data file for the spectrum corresponding to the redshift-dependent UVB, you need to specify the directory that contains all of the data files for the different redshifts, rather than one individual file. This directory will also need to contain an HDF5 file called ``redshifts.hdf5`` that gives the redshifts of the individual cross sections files (you can find an example of this in the ``HM12_cross_sections`` directory in the ``chimes-data`` repository). The cross sections files themselves need to be named as ``z%.3f_cross_sections.hdf5`` where the redshift is given to three decimal places (i.e. as ``%.3f``). CHIMES will then be able to find the cross sections files that bracket the current redshift. 

* ``use_redshift_dependent_eqm_tables`` - Set this parameter to ``1`` to interpolate the equilibrium abundance tables to the current redshift as well. This is important if you are using equilibrium tables that use the redshift-dependent UVB, since they will depend on redshift. If you are not using the equilibrium abundance tables, this parameter can be set to ``0``.

* ``EqAbundanceTablePath[500]`` - If you are using a redshift-dependent UVB and ``use_redshift_dependent_eqm_tables == 1``, this path needs to specify the directory containing the equilibrium table files for each redshift, rather than an individual file itself. These need to be given at the same redshifts as the corresponding cross sections data files. The equilibrium table files need to be named as ``z%.3f_eqm.hdf5``, where the redshift is given to three decimal places (i.e. as ``%.3f``).

* ``redshift`` - The current redshift. As this is defined in ``globalVariables``, rather than ``gasVariables``, it only needs to be set once at the beginning of each hydrodynamic time-step, rather than for each particle/cell individually.

* ``reionisation_redshift`` - The redshift of reionisation. At earlier redshifts, we set the UVB to zero. Also, if the two redshift bins that bracket the current redshift also bracket the redshift of reionisation, we do not interpolate. We just take either the earlier or the later table, if we are before or after reionisation, respectively. This avoids problems with interpolating between tables that are before and after reionisation. 

Interpolate UVB
^^^^^^^^^^^^^^^

At the beginning of each hydrodynamic timestep you will need to call the ``interpolate_redshift_dependent_UVB(*myGlobalVars)`` function, where ``*myGlobalVars`` is a pointer to the instance of the ``globalVariables`` structure. If this is the first time that this function has been called, it will load the cross sections tables that bracket the current redshift. Otherwise, it will compare the current redshift to the two tables that have already been loaded. If the current redshift is less than the lower redshift table, that table will be copied over to the high redshift table and the new lower redshift table will be read in, so that they again bracket the current redshift.

It will then interpolate the cross sections tables and, if necessary, the equilibrium abundance tables to the current redshift. If the current redshift is before reionisation, the strength of the UVB is set to zero. We therefore just set the tables the first time this routine is called, and then skip the interpolation on later calls until after reionisation, since the rates will just be zero. Also, as noted above, if the two tables bracket the redshift of reionisation, we do not interpolate them. We instead just take the earlier or the later tables, if we are before or after reionisation, respectively. 

References
^^^^^^^^^^

| `Faucher-Giguere (2020) <https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.1614F>`_ 
| `Haardt & Madau (2012) <https://ui.adsabs.harvard.edu/abs/2012ApJ...746..125H>`_ 
