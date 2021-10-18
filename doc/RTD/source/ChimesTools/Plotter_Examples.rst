.. Plotter Examples
   Alexander Richings, 16th March 2020

.. _PlotterExamples_label:

Examples
^^^^^^^^

The ``chimes-plotter/examples/`` directory in the ``chimes-tools`` repository contains example parameter files for running the ``chimes-plotter.py`` script. These can be used to plot the outputs from the corresponding Chimes Driver examples. 

plot_cooling_rates
""""""""""""""""""

This example is found in ``cooling_rates/plot_cooling_rates.param``. It plots the cooling and heating rates versus temperature for different densities and metallicities from the cooling tables. 

.. code-block:: none 

  plotter_mode  cooling_rates

Defines the mode in which Chimes Plotter will be run, i.e. to plot the cooling and heating rates from the cooling tables. 

.. code-block:: none 

  input_file   cooling_rates.hdf5
  output_file  cooling_plot.png

Name of the input HDF5 cooling table file and the output file for the plot. 

.. code-block:: none 

  log_y_min -30.0 
  log_y_max -20.0

Sets the minimum and maximum of the y-axis. The rates are plotted logarithmically in units of erg cm^3 s^-1. 

.. code-block:: none 

  x_variable  temperature
  log_x_min   2.0 
  log_x_max   9.0 

Plots temperature on the x-axis, from log10(T [K]) = 2.0 to 9.0. 

.. code-block:: none 

  line_variable  metallicity
  N_lines        3
  line_var_min   -2.0
  line_var_max   0.0

Different line styles represent metallicity. We plot for three metallicities, from log10(Z / Zsol) = -2.0 to 0.0. 

.. code-block:: none 

  panel_variable density
  N_panels  9
  panel_var_min -4.0
  panel_var_max 4.0

Different panels show different densities. We plot 9 panels, from log10(nH [cm^-3]) = -4.0 to 4.0. 

plot_eqm_abundances
"""""""""""""""""""

This example is found in ``eqm_table/plot_eqm_abundances.param``. It plots the equilibrium abundances versus temperature for different densities and metallicities from the equilibrium abundance tables. 

.. code-block:: none 

  plotter_mode  eqm_table

Defines the mode in which Chimes Plotter will be run, i.e. to plot the equilibrium abundances from the equilibrium abundance tables. 

.. code-block:: none 

  input_file   eqm_table.hdf5
  output_file  eqm_plot.png

Name of the input HDF5 equilibrium abundance table file and the output file for the plot. 

.. code-block:: none 

  species  elec
  species  HI
  species  HII 
  species  H2
  species  HeI 
  species  HeII 
  species  HeIII 

Defines which species to include in the plot. 

.. code-block:: none 

  log_y_min  -6.0 
  log_y_max  0.5

Sets the minimum and maximum of the y-axis. The abundances are plotted logarithmically relative to the corresponding element, i.e. log10(n_i / n_element). 

.. code-block:: none 

  x_variable  temperature
  log_x_min   2.0 
  log_x_max   9.0 

Plots temperature on the x-axis, from log10(T [K]) = 2.0 to 9.0. 

.. code-block:: none 

  line_variable  metallicity
  N_lines        3
  line_var_min   -2.0
  line_var_max   0.0

Different line styles represent metallicity. We plot for three metallicities, from log10(Z / Zsol) = -2.0 to 0.0. 

.. code-block:: none 

  panel_variable density
  N_panels  9
  panel_var_min -4.0
  panel_var_max 4.0

Different panels show different densities. We plot 9 panels, from log10(nH [cm^-3]) = -4.0 to 4.0. 

plot_species_evolution
""""""""""""""""""""""

This example is found in ``grid_noneq_evolution/plot_species_evolution.param``. It plots the evolution of the species abundances for different densities and metallicities. 

.. code-block:: none 

  plotter_mode  noneq_evolution_species 

Defines the mode in which Chimes Plotter will be run, i.e. to plot the abundance evolution. 

.. code-block:: none 

  input_file   grid_noneq_evolution.hdf5
  output_file  species_evolution_plot.png

Name of the input HDF5 file and the output file for the plot. 

.. code-block:: none 

  species  elec
  species  HI
  species  HII 
  species  H2
  species  HeI 
  species  HeII 
  species  HeIII 

Defines which species to include in the plot. 

.. code-block:: none 

  log_x_min  0.0
  log_x_max  3.0

Sets the minimum and maximum of the x-axis. The time is plotted logarithmically in Myr. 

.. code-block:: none 

  log_y_min  -6.0 
  log_y_max  0.5

Sets the minimum and maximum of the y-axis. The abundances are plotted logarithmically relative to hydrogen, i.e. log10(n_i / n_Htot). 

.. code-block:: none 

  fixed_variable   temperature 
  fixed_var_value  6.0 

All lines and panels will be plotted for the same initial temperature of log10(T [K]) = 6.0. 

.. code-block:: none 

  line_variable  metallicity
  N_lines        3
  line_var_min   -2.0
  line_var_max   0.0

Different line styles represent metallicity. We plot for three metallicities, from log10(Z / Zsol) = -2.0 to 0.0. 

.. code-block:: none 

  panel_variable density
  N_panels  9
  panel_var_min -4.0
  panel_var_max 4.0

Different panels show different densities. We plot 9 panels, from log10(nH [cm^-3]) = -4.0 to 4.0. 

plot_temperature_evolution
""""""""""""""""""""""""""

This example is found in ``grid_noneq_evolution/plot_temperature_evolution.param``. It plots the evolution of the temperature for different densities and metallicities. 

.. code-block:: none 

  plotter_mode  noneq_evolution_temperature 

Defines the mode in which Chimes Plotter will be run, i.e. to plot the temperature evolution. 

.. code-block:: none 

  input_file   grid_noneq_evolution.hdf5
  output_file  temperature_evolution_plot.png

Name of the input HDF5 file and the output file for the plot. 

.. code-block:: none 

  log_x_min  0.0
  log_x_max  3.0

Sets the minimum and maximum of the x-axis. The time is plotted logarithmically in Myr. 

.. code-block:: none 

  log_y_min  0.5
  log_y_max  8.5 

Sets the minimum and maximum of the y-axis. The temperature is plotted logarithmically in K. 

.. code-block:: none 

  N_T_init         7
  log_T_init_min   2.0
  log_T_init_max   8.0

Plot the evolution for 7 different initial temperatures, from log10(T_init [K]) = 2.0 to 8.0. Each initial temperature will be represented by a different line colour in each panel and for each line style. 

.. code-block:: none 

  line_variable  metallicity
  N_lines        3
  line_var_min   -2.0
  line_var_max   0.0

Different line styles represent metallicity. We plot for three metallicities, from log10(Z / Zsol) = -2.0 to 0.0. 

.. code-block:: none 

  panel_variable density
  N_panels  9
  panel_var_min -4.0
  panel_var_max 4.0

Different panels show different densities. We plot 9 panels, from log10(nH [cm^-3]) = -4.0 to 4.0. 





  






