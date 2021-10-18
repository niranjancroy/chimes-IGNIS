.. Parallelisation
   Alexander Richings, 23rd March 2020

.. _Parallelisation_label:

Parallelisation
---------------

When the CHIMES solver is called on a given particle/cell, it is run as a serial code, i.e. it is only using one CPU to integrate the chemistry and cooling for that one particle/cell. In a hydrodynamics code, the problem is typically parallelised by dividing the particles/cells between the available CPUs, using multi-threading and/or MPI. The chemistry and cooling integration itself depends only on local quantities, i.e. the properties of the particle/cell itself. This means that it does not require any communication between particles/cells, which in some respects makes the parallelisation fairly straightforward. 

The main complication is that the computational cost of the chemistry solver can be very different for different particles/cells. If the gas is already close to chemical and thermal equilibrium, the integration will be very quick. However, if it is far from equilibrium, the time-step will need to be sub-cycled many times to perform the integration, which will make it much more expensive. This is particularly problematic when using MPI to parallelise the simulation, because you would typically divide the particles/cells between the MPI tasks first and then do the chemistry and cooling integration on each one. Since it is difficult to know beforehand which particles/cells will be expensive, this means that a single MPI task can easily end up with most of the expensive particles/cells. This leads to poor work load balancing between MPI tasks, even if they all have a similar number of particles/cells.

I find that the best approach is to use a hybrid MPI + multi-threading method, where you would run 1 MPI task on each node and then each MPI task creates multiple threads equal to the number of CPUs per node. When you divide the loop over active particles between multiple threads (for example, with an OPENMP for loop), it can dynamically assign each particle/cell to the threads as it loops through them. In other words, a thread would integrate the chemistry and cooling for one particle, then when it is finished it would look for the next particle/cell in the loop. This means that, if one thread gets an expensive particle/cell, the other threads can continue with all of the other particles/cells on that node. This makes it more efficient than a pure MPI code, where you would typically need to divide up all of the particles/cells between the CPUs before you start the integration. 

MPI is still needed to run the simulation across multiple nodes. However, this hybrid approach reduces how many MPI tasks you need, which tends to improve the work load balancing and scaling.

Finally, if running with MPI, note that each MPI task will need to initialise CHIMES, so that they each have a copy of the reaction data tables etc. This can be done by having each task call ``init_chimes(*myGlobalVars)`` (see the :ref:`InitChimes_label` section). 
