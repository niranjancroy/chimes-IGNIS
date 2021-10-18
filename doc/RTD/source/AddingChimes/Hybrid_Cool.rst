.. Hybrid Cooling
   Alexander Richings, 26th March 2020

.. _HybridCool_label:

Hybrid Cooling
--------------

The radiative cooling and heating rates are typically computed from the non-equilibrium abundances of the ions and molecules in the CHIMES network. However, in some cases you may want to include additional cooling or heating channels that are not already in the network.

If the additional cooling or heating rate is constant over the course of the chemistry and cooling integration, you can use the ``gasVariables.constant_heating_rate`` variable (see the :ref:`DataStructures_label` section for details). For example, when you incorporate CHIMES into a hydrodynamics code there will be a cooling term due to adiabatic cooling (or heating) from the hydro solver. You could apply this in an operator-split fashion, where you update the internal energy of the particle/cell from the adiabatic cooling and then integrate the chemistry and radiative cooling. Alternatively, you could pass this adiabatic term to the ``constant_heating_rate`` variable, which would then add it on to the radiative cooling and heating rates, so that the adiabatic and radiative cooling is integrated together over the hydrodynamic time-step.

However, in some cases the additional cooling and heating rates may vary over the course of the integration. One example is if you only want to evolve some of the elements in non-equilibrium in the network, but still include the cooling from the other elements using pre-computed cooling tables in chemical equilibrium. Since these rates would vary with temperature, they would need to be updated throughout the integration.

To this end, we have provided an option for the User to implement their own additional cooling function that can be added within the CHIMES cooling routines. We call this "hybrid" cooling, which refers to using a hybrid approach of non-equilibrium for some elements and equilibrium for others - however, the User may have other applications for this option as well.

This hybrid cooling option is controlled through variables in the ``globalVariables`` and ``gasVariables`` structures. These are summarised in the :ref:`DataStructures_label` section, but we discuss these in more detail below.

Global Variables
^^^^^^^^^^^^^^^^

* ``hybrid_cooling_mode`` - An integer flag that switches on the Hybrid Cooling option. If this is set to ``0``, the hybrid cooling function (below) will not be used. If it is set to ``1``, then the hybrid cooling function will be called from within the CHIMES cooling routines and added on to the radiative cooling rates.

* ``*hybrid_cooling_fn`` - This is a function pointer that points to the User-defined hybrid cooling function. The idea is that the User can create their own function to define the additional cooling rates that they want to include, and then they set this function pointer to point to their function. CHIMES can then call their function from within the CHIMES cooling routines. This function must take two arguments, a pointer to a ``gasVariables`` structure and a pointer to a ``globalVariables`` structure, and return a double giving the net heating rate in units of ``erg cm^3 s^-1`` (i.e. the rate per unit volume divided by ``nH^2`` where ``nH`` is the total hydrogen number density), where heating is positive and cooling is negative. If the Hybrid Cooling option is not being used (i.e. if ``hybrid_cooling_mode == 0``), this pointer will be set to ``NULL`` by default.

* ``*hybrid_data`` - A void pointer that can be used to point to a User-defined structure that contains additional data that will be needed by the hybrid cooling function. Since this function takes the ``globalVariables`` as an argument, it will then have access to this hybrid data. It will need to cast this void pointer to a pointer to the User-defined structure. The hybrid data in ``globalVariables`` needs to be data that will be the same for all gas particles/cells (for example, any additional cooling tables). For data that is specific to a given gas particle/cell, there is another ``*hybrid_data`` pointer within ``gasVariables`` (which will point to a different User-defined structure; see below).

* ``*allocate_gas_hybrid_data_fn`` - A function pointer that points to a User-defined function to allocate memory to the ``gasVariables.hybrid_data`` structure (for example, if it contains dynamically allocated arrays). This is the hybrid data in ``gasVariables``, not ``globalVariables``. We need this function for the ``gasVariables`` hybrid data because we will need to allocate memory for it every time we create an instance of ``gasVariables``. This lets us do it all from within ``allocate_gas_abundances_memory()``, which will also call this function if ``hybrid_cooling_mode == 1``. Since the ``globalVariables`` instance is only set up once, at the beginning of the simulation, this can just be done from within the hydrodynamics code. The ``*allocate_gas_hybrid_data_fn`` function needs to take one argument, a pointer to the ``gasVariables`` structure, and return ``void``.

* ``*free_gas_hybrid_data_fn`` - A function pointer that points to a User-defined function to free memory from the ``gasVariables.hybrid_data`` structure. This is the hybrid data in ``gasVariables``, not ``globalVariables``. We need this function for the ``gasVariables`` hybrid data because we will need to free memory from it every time we destroy an instance of ``gasVariables``. This lets us do it all from within ``free_gas_abundances_memory()``, which will also call this function if ``hybrid_cooling_mode == 1``. The ``*free_gas_hybrid_data_fn`` function needs to take one argument, a pointer to the ``gasVariables`` structure, and return ``void``.

Gas Variables
^^^^^^^^^^^^^

* ``*hybrid_data`` - A void pointer that can be used to point to a User-defined structure that contains additional data that will be needed by the hybrid cooling function. This is different from the ``*hybrid_data`` structure in ``globalVariables``. As noted above, the ``*hybrid_data`` in ``gasVariables`` should be used for data that is specific to a given gas particle/cell. Since the hybrid cooling function takes the ``gasVariables`` as an argument, it will then have access to this hybrid data. It will need to cast this void pointer to a pointer to the User-defined structure. You will also need to provide functions to allocate and free memory for this structure (for example, if it contains dynamically allocated arrays), via ``*allocate_gas_hybrid_data_fn`` and ``*free_gas_hybrid_data_fn`` in ``globalVariables``. 
