.. Generate Cross Sections
   Alexander Richings, 3rd June 2020

.. _GenerateCrossSections_label: 

Generate Cross Sections
-----------------------

The ``generate_cross_sections.py`` python script can be used to generate the cross sections files needed for a given UV spectrum, which includes the averaged cross sections and photoheating rates given the shape of the UV spectrum (see the :ref:`ChimesCrossSections_label` section for details of what is included in these files). The ``chimes-data`` repository already contains cross sections files for several commonly used UV spectra, but if you want to use your own spectrum you will need to create a cross sections file from the spectrum using ``generate_cross_sections.py``.

To run this script, you will need to pass it a text file containing the various parameters, as discussed in the :ref:`CrossSectionsParameters_label` section below. The script can then be run on a single CPU as follows:

``python generate_cross_sections.py parameter_file``

It can also be run on multiple CPUs using MPI. For example, to run it on 4 CPUs you can use the following:

``mpirun -np 4 python generate_cross_sections.py parameter_file`` 

.. toctree::
    :maxdepth: 1

    Cross_Sections_Parameters
    Cross_Sections_Example

