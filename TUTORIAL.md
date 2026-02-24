# GPFPlot Tutorial

Version 0.2.1, February 2025

Before starting, make sure you have read the [README.md](README.md) file.

## 1. Simple usage of the program

1.1 Let us use the provided example file. In the examples directory, the files
`C2H4.inp`, `C2H4.out`, `C2H4.gpfep`, and `C2H4.dmat` are from a GAMESS / VB2000
calculation of the C2H4 molecule at GVB-PP/cc-pVTZ level. Move the files
`*.out`, `*.dmat` and `example1.gpfplot` to the parent directory (where
`gpfplot.py` is located).

1.2 Open the file example1.gpfplot and examine its content. All possible
GPFPlot variables are written in the file, although it is not
required (some are just defined with the default values), and some
lines are commented.

For example, `plane` is set to `xy`, which is coincidentally the plane of
the molecule, as it is defined in the GAMESS / VB2000 files.

The variable mode is set to `ORB 1`, that is, GPFPlot will plot the
first orbital of the GVB-PP function (in the case, it is one of the core
1s orbitals of the carbon atom).

1.3 Run `gpfplot.py` in the terminal. You shall see the initial text of the
program. Afterwards, it prints the main input files (`out_file` and
`dens_file`) and prints the selected mode of the program. The program will
then read the variables and calculate the orbital grids. Finally, the
program will plot the selected orbital. A new window should open. This
is the standard MatPlotLib plot window. There are some built-in options
from MatPlotLib available, as zooming the graph or save it as an image.
Notice, however, that part of the plot is empty, because it is off scale.

1.4 The values of `min_c` and `max_c` are not appropriate in this case. Changing
`auto_c` to `yes` should fix this problem, it will calculate the best
scale automatically. Close the plot window, edit the input file and run
GPFPlot again.

1.5 You can now change the mode variable in order to see the other orbitals
(`ORB 2`, `ORB 3`, etc.), or the HF orbitals from GAMESS (`HFORB 1`,
`HFORB 2`, etc.). Note that the orbitals of the pi bond (ORBs 13 and 14
or HFORBs 8 and 9) are not visible, because the xy plane is the nodal
plane of these orbitals.

1.6 You will be able to get a top view of those pi orbitals if set the offset
variable for, e.g. `0.30` Angstrom above the xy plane. When this is done,
the atom labels are not shown on the plot anymore (because the atoms are
not on the plotted plane), unless you set the `draw_atom` variable to `all`. 

1.7 Alternatively, you can change the plane to `xz` for better viewing the
pi orbitals.

1.8 Try to change the plot settings variables to see their effects. Increase
the value of `nconts` (I personally prefer about 30). Turn off `c_fill` and
uncomment `color_l`. Now the contour lines are not filled and are colored
themselves. Try changing draw_name to `no`. Try out other colormaps.
 
1.9 Finally, try to use the density program modes: set mode to `INT 3 4` in
order to visualize the interference of the C-C sigma bond. Try to
visualize other density partitions between other orbital pairs.

  
## 2. The `PROMPT` mode

The `PROMPT` mode is convenient for exploring a just finished calculation.
You can visualize various orbitals and densities without having to close
and rerun the program.

2.1 Change the mode to `PROMPT`. The program will initialize and, instead of
showing a plot window, it will wait for you type a command.

2.2 If you type, for example, `ORB 1` and press enter, it will plot the
orbital VB #1. You can type `ORB N`, `HFORB N`, `TOT X Y Z...`, etc., in
order to view any orbital or density with the geometric and plot
parameters previously selected.

2.3 You can type `exit` for closing PROMPT mode.

2.4 You can change the `# Plot Settings` variables defined in the input
file during the PROMPT session (`c_fill`, `color_f`, `c_lines`, `color_l`, 
`c_label`, `min_c`, `max_c`, `nconts`, `auto_c`, `p_unit`, `draw_atom`, 
and `draw_name`). Just type `variable=value`, for example, `auto_c=yes`.
You can set more than one variable at once, separated by blankspaces. For
example, `c_fill=no color_l=Spectral nconts=30 title="Alternative plot"`.

*Warning: if you choose `PROMPT` mode with a relatively high `gridp` number,*
*it will consume a lot of RAM, because the program calculates all*
*orbital grids and stores them in the memory.*

## 3. Saving Files

By defining the variables `save_txt` or `save_png` in the input file, it
becomes possible to save your plots as text files or PNG images.

3.1 Uncomment the two last lines in the example1.gpfplot file and change the
mode variable to, for example, `ORB 3`. The program will not only show
the plot window, but also generate two files.

The text file can be read by the script `gpfplot_figure.py`, but its
format can be easily readable to other plotting programs of your choice.
The first 10 lines of the file are comments specifying the grid
properties. The next lines contain the plotted orbital or density in
the form of a *N* x *N* points grid, where *N* is the gridp variable.


## 4. The `gpfplot_figure.py` script

If you want to make a figure with two or more plots, the `gpfplot_figure.py`
script may be helpful. Its structure is somewhat more raw than GPFPlot
itself and may have bugs.

This script is more like a template which must be manually customized
for every figure you want to create. You can make a copy of the original
script and move it to a different directory. In order to work in a 
different directory, you have to uncomment the lines 14-16 of the script
and edit the variable GPFPATH with the path of the GPFPlot installation.

4.1 Let us use the files `example2_fig1.gpfplot`, `example2_fig2.gpfplot`, and
`example2_fig3.gpfplot`. Move them from `examples` to its parent directory
and run them on GPFPlot. It will generate `*.txt` grid files for the
density partitioning of the sigma C-C bond of the C2H4 molecule
(recall that `C2H4.out` and `C2H4.dmat` must be on the same directory).

4.2 The `gpfplot_figure` script has no input files, it has to be manually
edited instead. Open it in a text editor and examine the section "read
variables" (lines 24-47).The following variables are relevant:

`figure_name`: the name of the final image generated by the script
(without extension).

`n_cols`: number of plots arranged horizontally. 

`n_rows`: number of plots arranged horizontally.

For example, if `n_cols` = 2 and `n_rows` = 2, the image generated will be
a grid of 2 x 2 = 4 plots (and 4 files must be given as input below).

`p_list`: List of `*.txt` files to be read and plotted. The program will
read from the `*.txt file` the locations of the `*.gpfplot` and `*.out` files
used initially, so they cannot be moved.

`titles`: The title of each plot. You can use a LaTeX-like language to
write Greek letters and mathematical expressions. For more information,
see [https://matplotlib.org/3.1.0/tutorials/text/mathtext.html].

`suptitle`: The title of the whole figure. Also admits math text options.

`p_size`: the print size of each plot in inches. I recommend do not change
this variable, because the font size of titles and labels could do not
match with other sizes.

`DPI`: resolution of the image in dots per inch.

`EPS`: option to generate a *.eps vector image.

4.3 The script is already set with the right variables to generate the plot.
If you run the script, it first will show a plot window in order to
check if everything is correct. Sometimes a little manual adjustment is
necessary (probably this time, as I intentionally messed up some
parameters for didactical purposes :) ). In the plot window, you can click
the button "configure subplots" in order to adjust the positions of the plot
box, axis titles, etc. Probably it will fit better by sliding "left" to 0.06,
"bottom" to 0.16 and "top" to 0.86. You can set these values in the "Manual
adjustments" section of the script (lines 137-142) and rerun the script
to obtain the correct image.
