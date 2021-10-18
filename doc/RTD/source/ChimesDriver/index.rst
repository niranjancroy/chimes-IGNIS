.. CHIMES Driver
   Alexander Richings, 3rd March 2020 

.. _ChimesDriver_label: 

CHIMES Driver
=============

This section describes the ``chimes-driver.py`` python script that can be used to run the stand-alone CHIMES module on a variety of applications, for example to compute equilibrium abundances from a simulation snapshot in post-processing, to follow the non-equilibrium chemical equilibrium on a grid of initial temperatures, densities and metallicities, or to create tables of cooling rates or equilibrium abundances. 

To run ``chimes-driver.py``, you need to pass it a parameter file as an argument. This is a text file defining the various parameters, which are described in the :ref:`ChimesDriverParam_label` section below. To run on a single CPU, you can use: 

``python chimes-driver.py parameter_file`` 

You can also run on multiple CPUs with MPI. For example, to run on 4 CPUs: 

``mpirun -np 4 python chimes-driver.py parameter_file``

By default, ``chimes-driver.py`` will assume that you built the CHIMES library in double precision by passing it the ``-DCHIMES_USE_DOUBLE_PRECISION`` compiler flag in the Makefile. If instead you built it in single precision, i.e. without the ``-DCHIMES_USE_DOUBLE_PRECISION`` compiler flag, you will need to run ``chimes-driver.py`` with an additional ``--single-precision`` argument, which must come after the parameter file, i.e.:

``mpirun -np 4 python chimes-driver.py parameter_file --single-precision``

You will also need to make sure that the Sundials library was built in single precision. 

.. toctree:: 
    :maxdepth: 1

    Parameters
    Inputs
    Outputs
    UV_Field
    Shielding
    Dust_Depletion
    Examples
