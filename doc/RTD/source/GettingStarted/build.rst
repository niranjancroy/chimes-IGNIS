.. Building CHIMES
   Alexander Richings, 19th February 2020 

.. _build_label:

Building CHIMES
--------------- 

The main ``chimes`` repository includes a template Makefile that can be used to build CHIMES as a stand-alone C library. 

First, you will need to create a local ``Makefile`` from the template using ``cp Makefile_template Makefile``. You may need to edit your local copy of the ``Makefile`` for your system. In particular, if it cannot automatically locate any system-installed HDF5 and Sundials libraries (for example, if you installed your own Sundials libraries), then you will need to specify the paths to these libraries via the ``SUNDIALS_PATH``, ``HDF5_INCL_PATH`` and ``HDF5_LIB_PATH`` variables at the beginning of the Makefile.

There are also compiler flags that you can pass to the compiler via the ``CFLAGS`` variable in the Makefile:

* ``-DCHIMES_USE_DOUBLE_PRECISION`` - Build CHIMES in double precision. If you do not include this flag, CHIMES will be built in single precision. Note that you will need to ensure that the Sundials library is built with the same precision as CHIMES. For new users, we recommend using double precision. Single precision can be faster in some circumstances, but it can also be slower. In general, when integrating the chemistry over a short period of time (for example, if you are following the non-equilibrium chemical evolution in small steps), single precision tends to perform better. However, when integrating for a long period of time (~1 Gyr or more), for example if you are just evolving the system to the final equilibrium state, single precision tends to be a lot slower than double precision, as it has to sub-cycle the integration with more steps (although they both reach the same solution). If you wish to explore using single precision, you will need to test it carefully to see whether single or double is fastest for your application.
  
* ``-DCHIMES_USE_GNU_SOURCE`` - This tells CHIMES to use the GNU extensions, which will allow it to use the ``exp10()`` function that is provided by the GNU extensions. If you do not include this compiler option, CHIMES will define its own ``exp10()`` function. If you are unsure whether your system supports GNU extensions, it is best to avoid this option. 

Once you have set up the Makefile for your system, you can then build CHIMES by running:

.. code-block:: bash 

  make 

This will build the CHIMES library, ``libchimes.so``, in the root ``chimes/`` directory. 
