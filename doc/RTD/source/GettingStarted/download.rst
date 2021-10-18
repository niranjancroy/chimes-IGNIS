.. Downloading CHIMES 
   Alexander Richings, 19th February 2020 

.. _Download_label: 

Downloading CHIMES
------------------

The CHIMES code and data files can be downloaded from the following Bit Bucket repositories using ``git``. These repositories are all publicly accessible, so you do not need to create an account to download them. 

CHIMES source code
^^^^^^^^^^^^^^^^^^

The CHIMES source code can be downloaded from the `main CHIMES repository <https://bitbucket.org/richings/chimes>`_. As it is a ``git`` repository, it can be cloned with the following command: 

.. code-block:: bash 

  git clone https://bitbucket.org/richings/chimes 

This contains the source code for the CHIMES module. It can be built as a stand-alone library, or the source code can be incorporated directly into your own hydrodynamics code. 

CHIMES data files
^^^^^^^^^^^^^^^^^

To run CHIMES, you will need the CHIMES data files that contain the various reaction rate coefficients, photoionisation cross sections etc. These can be downloaded from the `chimes-data <https://bitbucket.org/richings/chimes-data>`_ repository. This repository uses Git LFS (Large File Storage) to store the HDF5 data files, so you will need to install ``git-lfs`` on your system before you can clone it (see `this page <https://github.com/git-lfs/git-lfs/wiki/Installation>`_ for details). Once you have installed ``git-lfs``, the ``chimes-data`` repository can then be cloned with ``git`` using the following command: 

.. code-block:: bash 

  git clone https://bitbucket.org/richings/chimes-data

You will need to make a note of the path to your local copy of chimes-data, as you will need to pass this as a parameter to the CHIMES module. 

CHIMES Driver
^^^^^^^^^^^^^

If you build CHIMES as a stand-alone C library, it can be used from within a Python script using the ``ctypes`` python module. We have put together such a python script, ``chimes-driver.py``, to demonstrate how to do this. CHIMES Driver can be used to apply the stand-alone CHIMES module to a variety of applications, including to run it on simulation snapshots in post-processing (either to compute the equilibrium chemical abundances of the gas particles in a snapshot, or to follow the non-equilibrium chemical evolution of the gas particles at fixed density); to produce equilibrium abundance tables and cooling tables as a function of density, temperature and metallicity; or to follow the non-equilibrium chemical evolution over a grid of physical parameters. 

CHIMES Driver can be downloaded from the `chimes-driver <https://bitbucket.org/richings/chimes-driver>`_ repository. This repository can be cloned with ``git`` using the following command: 

.. code-block:: bash 

  git clone https://bitbucket.org/richings/chimes-driver 

CHIMES Tools
^^^^^^^^^^^^

CHIMES Tools contains various tools and scripts that can be useful, such as plotting scripts and routines for parsing the output CHIMES abundance arrays. CHIMES Tools is not required for running the CHIMES module. It can be downloaded from the `chimes-tools <https://bitbucket.org/richings/chimes-tools>`_ repository, which can be cloned with ``git`` using the following command: 

.. code-block:: bash 

  git clone https://bitbucket.org/richings/chimes-tools 
