.. Plotter Parameters
   Alexander Richings, 11th March 2020

.. _PlotterParameters_label:

Parameters
^^^^^^^^^^

The following parameters can be specified in the parameter file passed to ``chimes-plotter.py``. If a parameter is not given in the parameter file, it uses a default value as defined in ``chimes-plotter/plotter_config.py``. 

+-------------------------------------+------------------------------------------------------------------------------+
| Parameter                           | Description                                                                  |
+=====================================+==============================================================================+
| ``plotter_mode``                    | | Defines what type of Chimes Driver output we are plotting from. Possible   |
|                                     | | options are:                                                               |
|                                     | | ``eqm_table`` - Plot the equilibrium abundance tables produced by the      |
|                                     | |     corresponding ``eqm_table`` driver mode from Chimes Driver.            |
|                                     | | ``cooling_rates`` - Plot the cooling and heating rates from the cooling    |
|                                     | |     rate tables produced by the ``cooling_rates`` driver mode from Chimes  |
|                                     | |     Driver.                                                                |
|                                     | | ``noneq_evolution_temperature`` - Plot the temperature evolution produced  |
|                                     | |     by the ``noneq_evolution`` driver mode from Chimes Driver.             |
|                                     | | ``noneq_evolution_species`` - Plot the evolution of the species abundances |
|                                     | |     produced by the ``noneq_evolution`` driver mode from Chimes Driver.    |
|                                     | | Note that all of the above can only be used with Chimes Driver outputs     |
|                                     | | that use ``IO_mode == grid``. If ``IO_mode == snapshot``, the layout of    |
|                                     | | of the data within the snapshot will depend on the type of simulation      |
|                                     | | that produced the snapshot (e.g. Gizmo, Arepo etc.), so we have not        |
|                                     | | included general plotting routines for snapshot outputs.                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``input_file``                      | | Name of the input HDF5 file that we will be plotting from                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``output_file``                     | | Name of the output file for the plots that will be produced.               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``species <name>``                  | | Name(s) of the species that will be included in the plot, where            |
|                                     | | ``<name>`` is a valid species name as defined in CHIMES (see the           |
|                                     | | :ref:`ReadChimes_label` section for instructions on how to get a list of   |
|                                     | | valid species names from ``read_chimes.py``). Multiple species can be      |
|                                     | | defined, with one species per line in the parameter file, e.g.:            |
|                                     | | ``species <name1>``                                                        |
|                                     | | ``species <name2>``                                                        |
|                                     | | This is only used if ``plotter_mode`` is set to ``eqm_table`` or           |
|                                     | | ``noneq_evolution_species``.                                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``x_variable``                      | | The variable to be plotted on the x-axis. Possible options are:            |
|                                     | | ``temperature``, ``density`` or ``metallicity``. This is only used in if   |
|                                     | | ``plotter_mode`` is set to ``eqm_table`` or ``cooling_rates``. For the     |
|                                     | | ``noneq_evolution_temperature`` and ``noneq_evolution_species`` modes,     |
|                                     | | time is plotted on the x-axis.                                             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_x_min``                       | | The minimum value on the x-axis. This is expressed as log10(T [K]) for     |
|                                     | | temperature, log10(nH [cm^-3]) for density, log10(Z / Zsol) for            |
|                                     | | metallicity (where Zsol = 0.0129 is solar metallicity), and                |
|                                     | | log10(time [Myr]) for time.                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_x_max``                       | | The maximum value on the x-axis. This is expressed as log10(T [K]) for     |
|                                     | | temperature, log10(nH [cm^-3]) for density, log10(Z / Zsol) for            |
|                                     | | metallicity (where Zsol = 0.0129 is solar metallicity), and                |
|                                     | | log10(time [Myr]) for time.                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_y_min``                       | | The minimum value on the y-axis. This is expressed as log10(n_i / n_elem)  |
|                                     | | for ``plotter_mode == eqm_table``, log10(Lambda [erg cm^3 s^-1]) for       |
|                                     | | ``plotter_mode == cooling_rates``, log10(n_i / n_Htot) for                 |
|                                     | | ``plotter_mode == noneq_evolution_species``, and log10(T [K]) for          |
|                                     | | ``plotter_mode == noneq_evolution_temperature``.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_y_max``                       | | The maximum value on the y-axis. This is expressed as log10(n_i / n_elem)  |
|                                     | | for ``plotter_mode == eqm_table``, log10(Lambda [erg cm^3 s^-1]) for       |
|                                     | | ``plotter_mode == cooling_rates``, log10(n_i / n_Htot) for                 |
|                                     | | ``plotter_mode == noneq_evolution_species``, and log10(T [K]) for          |
|                                     | | ``plotter_mode == noneq_evolution_temperature``.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``panel_variable``                  | | The variable that is varied from one panel to another. Possible options    |
|                                     | | are: ``temperature``, ``density`` or ``metallicity``. If                   |
|                                     | | ``plotter_mode == noneq_evolution_temperature``, you cannot select         |
|                                     | | ``temperature`` here, because different initial temperatures are           |
|                                     | | represented by different line colours.                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``N_panels``                        | | The number of panels to include in the plot.                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``panel_var_min``                   | | The minimum value of the variable that is varied from one panel to         |
|                                     | | another, which will be used for the top-left panel. This is expressed as   |
|                                     | | log10(T [K]) for temperature, log10(nH [cm^-3]) for density and            |
|                                     | | log10(Z / Zsol) for metallicity (where Zsol = 0.0129 is solar              |
|                                     | | metallicity).                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``panel_var_max``                   | | The maximum value of the variable that is varied from one panel to         |
|                                     | | another, which will be used for the bottom-right panel. This is expressed  |
|                                     | | as log10(T [K]) for temperature, log10(nH [cm^-3]) for density and         |
|                                     | | log10(Z / Zsol) for metallicity (where Zsol = 0.0129 is solar              |
|                                     | | metallicity).                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``line_variable``                   | | The variable that is represented by different line styles (solid, dashed   |
|                                     | | etc). Possible options are: ``temperature``, ``density`` or                |
|                                     | | ``metallicity``. If ``plotter_mode == noneq_evolution_temperature``,       |
|                                     | | you cannot select ``temperature`` here, because different initial          |
|                                     | | temperatures are represented by different line colours.                    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``N_lines``                         | | The number of different line styles to include in the plot.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``line_var_min``                    | | The minimum value of the variable that is represented by different line    |
|                                     | | styles. This is expressed as log10(T [K]) for temperature,                 |
|                                     | | log10(nH [cm^-3]) for density and log10(Z / Zsol) for metallicity (where   |
|                                     | | Zsol = 0.0129 is solar metallicity).                                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``line_var_max``                    | | The maximum value of the variable that is represented by different line    |
|                                     | | styles. This is expressed as log10(T [K]) for temperature,                 |
|                                     | | log10(nH [cm^-3]) for density and log10(Z / Zsol) for metallicity (where   |
|                                     | | Zsol = 0.0129 is solar metallicity).                                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``fixed_variable``                  | | The variable that is held fixed for all panels and line-styles etc.        |
|                                     | | Possible options are: ``temperature``, ``density`` or ``metallicity``.     |
|                                     | | This is only used if ``plotter_mode == noneq_evolution_species``. This is  |
|                                     | | needed for this plotter mode because we cannot show variations of all      |
|                                     | | three variables with combinations of different panels and line styles      |
|                                     | | etc., so one of them has to be held fixed.                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``fixed_var_value``                 | | The value of the variable that is held fixed. This is expressed as         |
|                                     | | log10(T [K]) for temperature, log10(nH [cm^-3]) for density and            |
|                                     | | log10(Z / Zsol) for metallicity (where Zsol = 0.0129 is solar              |
|                                     | | metallicity).                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``N_T_init``                        | | Number of initial temperatures to plot evolution curves for in each panel  |
|                                     | | and for each line style. Different initial temperatures are represented    | 
|                                     | | by different line colours. Only used if                                    |
|                                     | | ``plotter_mode == noneq_evolution_temperature``.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_T_init_min``                  | | The minimum initial temperature to be plotted, expressed as                |
|                                     | | log10(T_init [K]). Only used if                                            |
|                                     | | ``plotter_mode == noneq_evolution_temperature``.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``log_T_init_max``                  | | The maximum initial temperature to be plotted, expressed as                |
|                                     | | log10(T_init [K]). Only used if                                            |
|                                     | | ``plotter_mode == noneq_evolution_temperature``.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

The following parameters control the appearance of the plots. They can generally be left at their default values, but they give the User the option to adjust their plots if required. 

+-------------------------------------+------------------------------------------------------------------------------+
| Parameter                           | Description                                                                  |
+=====================================+==============================================================================+
| ``figure_width``                    | | The width of the figure, in inches. Default: ``8.0``.                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``figure_height``                   | | The height of the figure, in inches. Default: ``5.0``.                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``subplot_left``                    | | Position of the left edge of the grid of subplots. Default: ``0.08``.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``subplot_right``                   | | Position of the right edge of the grid of subplots. Default: ``0.8``.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``subplot_bottom``                  | | Position of the bottom edge of the grid of subplots. Default: ``0.1``.     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``subplot_top``                     | | Position of the top edge of the grid of subplots. Default: ``0.95``.       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``line_width``                      | | Width of the plot lines. Default: ``1.8``.                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``tick_label_fontsize``             | | Font size for tick labels. Default: ``14``.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``axis_label_fontsize``             | | Font size for axis labels. Default: ``14``.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``panel_label_fontsize``            | | Font size for panel labels. Default: ``9``.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``legend_fontsize``                 | | Font size for legends. Default: ``9``.                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``panel_label_x``                   | | Relative horizontal position of the panel labels. Default: ``0.98``.       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``panel_label_y``                   | | Relative vertical position of the panel labels. Default: ``0.02``.         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``max_ticks``                       | | Maximum number of major ticks on each axis. Default: ``6``.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
