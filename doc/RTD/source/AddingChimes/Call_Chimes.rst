.. Call Chimes
   Alexander Richings, 23rd March 2020

.. _CallChimes_label:

Calling the CHIMES Solver
-------------------------

The CHIMES solver replaces the standard cooling routines in the hydrodynamics code. Every hydrodynamic time-step, at the point where you would do the standard cooling, you will need to do the following:

* Loop through all active gas particles/cells.

* Update the ``gasVariables`` structure with the particle/cell's current hydrodynamic quantities (density, temperature etc.). The ``gasVariables.hydro_timestep`` should be set to the length of the particle/cell's hydrodynamic time-step, in seconds. This tells CHIMES how long it needs to integrate the chemistry and cooling for.

* Call the ``chimes_network(*myGasVars, *myGlobalVars)`` function, where ``*myGasVars`` and ``*myGlobalVars`` are pointers to the ``gasVariables`` structure for the given particle/cell and the ``globalVariables`` structure, respectively. This is the function that will actually integrate the chemical abundances and, if ``gasVariables.ThermEvolOn == 1``, the temperature.

* If cooling was included in the above integration, you will then need to update the particle/cell's thermodynamic variable(s) in the hydrodynamics code (e.g. the internal energy or the entropy, depending on what is used in the hydro solver), using the final temperature that the particle/cell reached at the end of the time-step, as given in ``gasVariables.temperature``.

  
