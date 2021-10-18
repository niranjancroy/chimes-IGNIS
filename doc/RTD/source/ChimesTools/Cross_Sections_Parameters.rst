.. Cross Sections Parameters
   Alexander Richings, 3rd June 2020

.. _CrossSectionsParameters_label: 

Parameters
----------

The following parameters can be specified in the parameter file passed to ``generate_cross_sections.py``. If a parameter is not given in the parameter file, it uses a default value as defined in the ``read_parameter_file()`` function within ``generate_cross_sections.py``. 

+-------------------------------------+------------------------------------------------------------------------------+
| Parameter                           | Description                                                                  |
+=====================================+==============================================================================+
| ``chimes_main_data_path``           | | Path to the ``chimes_main_data.hdf5`` data file from the ``chimes-data``   |
|                                     | | repository.                                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``cross_sections_data_path``        | | Path to the directory containing additional data files needed by           |
|                                     | | ``generate_cross_sections.py`` (for example, the frequency-dependent       |
|                                     | | cross sections for each species). These can be found in                    |
|                                     | | ``chimes-tools/generate_cross_sections/data/``.                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``spectrum_file``                   | | Name of the text file that defines the shape of the UV spectrum. This      |
|                                     | | needs to contain a list of data points specifying log10(E [Ryd]) in the    |
|                                     | | first column and log10(J_nu [erg/s/sr/cm^2/Hz]) in the second column,      |
|                                     | | where E and J_nu are the photon energy and specific intensity              |
|                                     | | respectively. Lines that start with a ``#`` character are treated as       |
|                                     | | comments and are ignored. We include an example that can be used as a      |
|                                     | | template (see the :ref:`CrossSectionsExample_label` section).              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``output_file``                     | | Name of the output cross sections HDF5 file.                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

