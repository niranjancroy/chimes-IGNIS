.. Misc
   Alexander Richings, 26th March 2020

.. _Misc_label:

Miscellaneous
-------------
   
Advection and Diffusion
^^^^^^^^^^^^^^^^^^^^^^^

If you are using a particle-based hydro solver (e.g. SPH), the individual ions and molecules would typically just move around with the particles, so that there is no transfer of ions and molecules between particles. However, if you are using a grid-based code, gas will be advected between cells, and so individual ions and molecules will also need to be advected. Some particle-based codes have also implemented subgrid models for the diffusion of metals between particles (e.g. Hopkins 2017), in which case the ions and molecules can also diffuse between particles. In both cases, the advection or diffusion of individual ions and molecules can be treated in the same way as for the metals.

CHIMES Species Names
^^^^^^^^^^^^^^^^^^^^

The ``chimes_species_names`` array, which is declared in ``chimes_proto.h`` and defined in ``init_chimes.c``, gives an array of strings containing the names of all species in the full CHIMES network. This can be useful if you want to write out the species names as meta-data in your simulation outputs. If you are using a reduced network (i.e. if you have switched off individual elements), you can select the names of only the species included in the reduced network using the ``globalVariables.speciesIndices`` array (see the :ref:`DataStructures_label` section for details).

Doxygen Documentation
^^^^^^^^^^^^^^^^^^^^^

More detailed information on the functions and data structures in the CHIMES code can be found in the Doxygen documentation. To build the Doxygen documentation, you will need to ``cd`` into the ``doc/Doxygen`` directory in the ``chimes`` repository. From there, run ``doxygen`` (you will need to have ``doxygen`` installed on your system). This will build the documentation in ``doc/Doxygen/html/``. Open the ``index.html`` file in this sub-directory to open the main page of the documentation. 


References
^^^^^^^^^^

| `Hopkins (2017) <https://ui.adsabs.harvard.edu/abs/2017MNRAS.466.3387H>`_ 
