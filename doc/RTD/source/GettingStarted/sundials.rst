.. CHIMES Installing Sundials
   Alexander Richings, 19th February 2020 

.. _sundials_label: 

Installing Sundials
-------------------

CHIMES requires the CVODE library of Ordinary Differential Equation (ODE) solvers, from the Sundials package (v4.0 or later). If Sundials is not already installed on your system, it can be downloaded `here <https://computing.llnl.gov/projects/sundials/sundials-software>`_. 

Once you have downloaded Sundials from the above website, you can then untar it and build it. You will need ``cmake`` to build Sundials. The INSTALL_GUIDE.pdf included with the Sundials download describes this process in detail, but in brief these are the steps that you will need: 

.. code-block:: bash

  tar -zxvf sundials-5.1.0.tar.gz 
  cd sundials-5.1.0 
  mkdir build 
  cd build 
  cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir/ -DBUILD_ARKODE=OFF -DBUILD_CVODE=ON -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=ON -DCMAKE_C_FLAGS="-O2" -DEXAMPLES_ENABLE_C=OFF -DSUNDIALS_PRECISION=double ../
  make
  make install

In the above example, after untar'ing the sundials download we create a new ``build`` directory inside ``sundials-5.1.0``, and then ``cd build`` before running ``cmake``. This is needed because the build directory has to be different from the source directory. 

Also, we only need the CVODE library, so in this example we have switched off building the various other libraries that are included in the Sundials package, as these are not required by CHIMES.

You will need to ensure that both the Sundials library and CHIMES are built with the same precision (either single or double; see the :ref:`build_label` section). We recommend that new users build CHIMES in double precision, so in this example we use the ``-DSUNDIALS_PRECISION=double`` option to build Sundials in double precision as well. For single precision, you can use ``-DSUNDIALS_PRECISION=single``. Note that if Sundials is already installed on your system it is likely that it uses double precision, which is the default for Sundials. 

Finally, you will need to set the ``-DCMAKE_INSTALL_PREFIX`` path to a directory where you have write access. This is where the libraries will be installed. 

Once you have installed Sundials, you will need to add ``/path/to/install/dir/lib`` (or possibly ``/path/to/install/dir/lib64``, depending on your system) to your ``LD_LIBRARY_PATH`` environment variable, both when you compile CHIMES and when you run it. The command to set this variable depends on which shell you are using (to determine which shell you are using, you can run ``echo $SHELL`` in your terminal). The commands for two commonly used shells are as follows: 

For ``bash``, use: 

.. code-block:: bash

  export LD_LIBRARY_PATH=/path/to/install/dir/lib:$LD_LIBRARY_PATH 

For ``tcsh``, use: 

.. code-block:: tcsh

  setenv LD_LIBRARY_PATH /path/to/install/dir/lib:$LD_LIBRARY_PATH 

The above examples assume that ``LD_LIBRARY_PATH`` already exists, and prepends the Sundials path to it (which preserves any other paths that are already defined there). If it does not already exist, you will need to omit the ``:$LD_LIBRARY_PATH`` from the end of the command. 

You can run this command every time you compile or run CHIMES. Or, for convenience, you can simply add this command to your ``~/.bashrc`` or ``~/.tcshrc`` file (for your given shell), which will then define ``LD_LIBRARY_PATH`` every time you open a terminal. 
