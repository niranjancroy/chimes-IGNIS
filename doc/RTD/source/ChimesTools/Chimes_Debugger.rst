.. Chimes Debugger
   Alexander Richings, 16th March 2020

.. _ChimesDebugger_label: 

Chimes Debugger
---------------

``chimes-debugger`` is a simple C wrapper for the CHIMES module. It will set up and run the CHIMES chemistry solver on a single set of physical parameters (i.e. density, temperature, radiation field etc.). This is similar to what the ``chimes-driver.py`` python script does, except that it is much simpler as it only sets up the chemistry solver for a single set of parameters, which are hard-coded and are not read in from a parameter file. Also, ``chimes-debugger`` does not produce any outputs, it just runs through the chemistry solver and then exits. 

This was originally intended to make it easier to run CHIMES through a debugger such as GDB or DDT, rather than trying to debug a python code that imports CHIMES as a C-library. If you found an example where CHIMES crashes in ``chimes-driver.py``, you could then set up that example in ``chimes-debugger`` and run it through GDB/DDT directly. 

I have included ``chimes-debugger`` here in this repository as it is also useful as a simple example for how to incorporate CHIMES into a C code. This can be used as a template when adding CHIMES to your own hydro codes, for example. 

It is currently set up to run with a fixed temperature of 1e4 K, a density of 1 cm^-3 and solar metallicity, with the Black (1987) interstellar radiation field in the local solar neighbourhood and a shielding length set equal to the Jeans length. To change this set up, you will need to manually edit the parameters in ``chimes-debugger/main.c``. 

You will then need to create your own ``Makefile`` from the ``Makefile_template``. This is similar to the one from the ``chimes`` repository, and you may need to specify the paths to the Sundials and HDF5 libraries, depending on your system (see the :ref:`GettingStarted_label` section for details). You will also need to specify the path to your local copy of the main ``chimes`` repository, via the ``CHIMES_PATH`` variable in the ``Makefile_template``. 

You can then compile by running ``make`` as usual. This will create the ``chimes-debugger`` executable. You can then run it as: 

``chimes-debugger /path/to/chimes-data`` 

Where you need to provide the path to your local copy of the ``chimes-data`` repository as a command line argument. 

References
^^^^^^^^^^

| `Black (1987) <https://ui.adsabs.harvard.edu/abs/1987ASSL..134..731B>`_ 
