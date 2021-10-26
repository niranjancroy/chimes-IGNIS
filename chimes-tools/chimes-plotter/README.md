# chimes-plotter

This python script can be used to plot the outputs from the chimes-driver script. It is controlled by a parameter file (examples can found in the example_parameter_files directory). It is then run as follows:

python chimes-plotter.py parameters.param

The parameters are described below. Many of these can be omitted from the parameter file, in which case they will use the default values specified in plotter_config.py.

## Parameters
* ``plotter_mode`` -- This determines what type of output from chimes-driver we wish to plot, and should match the corresponding driver_mode from chimes-driver. Currently, we only support the eqm_table mode, which will plot the equilibrium abundances from an equilibrium table produced by chimes-driver.

* ``input_file`` -- Filename of the input data file, from chimes-driver.

* ``output_file`` -- Filename of the output plot.

* ``species XX`` -- This is used to specify which species to include in the plot, where XX gives the name of the species (e.g. HI). Multiple species can be specified, with one per line in the parameter file. Note that, currently, this only works with the full CHIMES network. If some elements have been switched off, we do not currently have the functionality to map from the species strings to the index in the reduced abundance array - this is something that we need to add to read_chimes.py in the near future.

* ``x_variable`` -- Specifies which variable of the table grid to plot on the x-axis. Must be one of: temperature, density, or metallicity.

* ``line_variable`` -- Specifies which variable of the table grid to represent with different line styles (i.e. solid, dashed etc.). Must be one of: temperature, density, or metallicity.

* ``N_lines`` -- Specifies how many different line styles to plot.

* ``line_var_min`` -- Specifies the minimum value of the line_variable to plot. Note that this is taken as a log.

* ``line_var_max`` -- Specifies the maximum value of the lin_variable to plot. Note that this is taken as a log. Then the spacing in line_variable is determined automatically from N_lines. 

* ``panel_variable`` -- Specifies which variable of the table grid to plot in different panels. Must be one of: temperature, density, or metallicity.

* ``N_panels`` -- Specifies the number of panels to include in the plot. 

* `` panel_var_min`` -- Specifies the minimum value of the panel_variable to plot. Note that this is taken as a log.

* ``panel_var_max`` -- Specifies the maximum value of the panel_variable to plot. Note that this is taken as a log. Then the spacing in panel_variable is determined automatically from N_panels.

* ``log_x_min`` -- Specifies the minimum of the x-axis, as a log.

* ``log_x_max`` -- Specifies the maximum of the x-axis, as a log.

* ``log_y_min`` -- Specifies the minimum of the y-axis, as a log.

* ``log_y_max`` -- Specifies the maximum of the y-axis, as a log.

### Appearance
The following parameters control the appearance of the figure. In general, they can be left to their default values, but if necessary the user can change them via the parameter file.

* ``figure_width`` -- Overall width of the figure, in inches. Default: 8.0. 

* ``figure_height`` -- Overall height of the figure, in inches. Default: 5.0. 

* ``subplot_left`` -- Relative position of the left-most edge of the subplots. Default: 0.08.

* ``subplot_right`` -- Relative position of the right-most edge of the subplots. Default: 0.8.

* ``subplot_bottom`` -- Relative position of the bottom-most edge of the subplots. Default: 0.1.

* ``subplot_top`` -- Relative position of the top-most edge of the subplots. Default: 0.95.

* ``line_width`` -- Width of the plot lines. Default: 1.8.

* ``tick_label_fontsize`` -- Font size used for the tick labels. Default: 14.

* ``axis_label_fontsize`` -- Font size used for the axis labels. Default: 14.

* ``panel_label_fontsize`` -- Font size used for the panel labels. Default: 9.

* ``legend_fontsize`` -- Font size used for the legends. Default: 9.

* ``panel_label_x`` -- Position of the right-most edge of each panel label, relative to the given panel. Default: 0.98.

* ``panel_label_y`` -- Position of the bottom-most edge of each panel label, relative to the given panel. Default: 0.02.

* ``max_ticks`` -- The maximum number of major ticks on each axis. Default: 6. 


