# GPFPlot
Version 0.2, October 2021

Created by David Sousa

Inspired by G. N. Freitas' DENSPLOT (non-published)

Also inspired by Brett Bode's [wxMacMolPlt](brettbode.github.io/wxmacmolplt)

Powered by Python / NumPy / MatPlotLib


## 1. WHAT IS GPF-PLOT?
--------------------
GPF-PLOT is a Python script which generates 2D plots of orbitals, electronic densities and GPF density partitioning associated to VB2000 output files.


## 2. CHANGELOG
------------
### Version 0.2, October 2021:
* Minor changes for adjusting to VB2000 3.0 / GAMESS 2021 R2 versions.
* Fixed error in normalization formula.
* Fixed parse_overlaps when there is more than 10 groups.
* Fixed parse_gpf when SPHER keyword is on.
* Added feature to save EPS file and customize Figure DPI.
### Version 0.1, July 2020:
* Initial release.


## 3. HOW TO RUN GPF-PLOT
----------------------
GPF-PLOT uses Python 3.x, so it needs to be installed in your machine. Additional requirements are the [NumPy](https://numpy.org) and [MatPlotLib](https://matplotlib.org) libraries. The current version of GPF-PLOT (0.2) was tested with Python 3.8.2, NumPy 1.17.4, and MatPlotLib 3.1.2.

If the above dependencies are installed, you should be able to run GPF-PLOT by entering the following command in your terminal (you must be in the directory which contains the script):

`python3 gpfplot.py INPUTFILE`

where `INPUTFILE` is a `*.gpfplot` valid input file.


## 4. PROGRAM MODES
----------------
Currently, GPF-PLOT has 6 program modes: `ORB`, `HFORB`, `QC`, `INT`, `TOT`, and `PROMPT`.

`ORB` and `HFORB` are respectively for plotting VB2000 and GAMESS output orbitals.

`QC`, `INT`, and `TOT` are for plotting the total electronic density (`TOT`) or its quasi-classical (`QC`) or interference (`INT`) densities, for a given set of orbitals.

`PROMPT` is an interactive prompt mode, in which you can ask for the program to plot orbitals and/or densities.

For more instructions about the modes, please see the [TUTORIAL.md](TUTORIAL.md) file.


## 5. HOW TO WRITE THE INPUTFILE
-----------------------------
Although `INPUTFILE` must have the `*.gpfplot` extension, it is just a plain text in which INPUT VARIABLES are set.

In order to set variables, type the name of the variable, followed immediately by an equal sign `=`, followed immediately by the variable value. Do not use spaces before or after the equal sign. Each line must contain only 1 variable (or none). Variables can be defined in any order. If a variable is defined more than once in the `INPUTFILE`, the
program reads ONLY the FIRST definition and ignore the rest. Lines starting with `#` are comments (ignored by the program). 


### 5.1 INPUT VARIABLES

#### 5.1.1 Input Files and program mode

VARIABLE  | DEFAULT VALUE |  DESCRIPTION                                
------------------------------------------------------------------------
`out_file` |       -       | Name of GAMESS/VB2000 output file with      
           |               | extension.                                  
------------------------------------------------------------------------
`dens_file`|       -       | Name of the density matrix file generated   
           |               | by VB2000.                                  
------------------------------------------------------------------------
`mode`     |   `PROMPT`    | Action which the program will perform.      
           |               |                                             
           |               | Possible values:                            
           |               |                                             
           |               |`PROMPT`: interactive prompt.                 
           |               |`ORB N`: plot Nth VB orbital.                 
           |               |`HFORB N`: plot Nth HF orbital.               
           |               |`QC X Y Z`:  plot quasi-classical density of 
           |               |             orbitals X, Y, Z, ...
           |               |`INT X Y Z`: plot interference density of    
           |               |             orbitals X, Y, Z, ...  
           |               |             (minimum of 2 orbitals).                          
           |               |`TOT X Y Z`: plot total density of orbitals
           |               |             X, Y, Z, ...
 -----------------------------------------------------------------------


#### 5.1.2 Geometry Settings

VARIABLE  | DEFAULT VALUE |  DESCRIPTION                                
------------------------------------------------------------------------
`plane`   |     `xy`      | Spatial plane where densities will be       
          |               | evaluated (options: xy, xz, yz, yx, zx, zy).
------------------------------------------------------------------------
`gridp`   |     `50`      | Number of points of the generated grid in   
          |               | each dimension (2D plot).                   
------------------------------------------------------------------------
`xlim`    |    `-1,1`     | Range for horizontal axis of the plot (Å).  
------------------------------------------------------------------------
`ylim`    |    `-1,1`     | Range for vertical axis of the plot (Å).    
------------------------------------------------------------------------
`offset`  |     `0.0`     | Value of the coordinate parallel to         
          |               | cartesian plane chosen for plotting.        
------------------------------------------------------------------------


#### 5.1.3 Plot Settings

VARIABLE  | DEFAULT VALUE |  DESCRIPTION                                
------------------------------------------------------------------------
`c_lines` |     `yes`     | Option for plotting isovalue lines (yes/no).         
------------------------------------------------------------------------
`c_fill`  |     `yes`     | Apply filled contour plot (yes/no).         
------------------------------------------------------------------------
`c_label` |     `no`      | Include labels int the isovalue lines       
          |               | (yes/no).                                   
------------------------------------------------------------------------
`color_f` |  `coolwarm`   | Colormap for filled contour plot. See       
          | (mode = ORB,  | section 6.1 for available colormaps.        
          |   HFORB or    |                                             
          |    PROMPT)    |                                             
          |  `seismic_r`  |                                             
          |  (mode = QC,  |                                             
          |  INT or TOT)  |                                             
------------------------------------------------------------------------
`color_l` |       -       | Colormap for isovalue lines. If not set,    
          |               | lines will be plotted in black. See section 
          |               | 6.1 for available colormaps.                
------------------------------------------------------------------------
`nconts`  |      `10`     | Number of contour isolines and/or color
          |               | changes in filled contour plot.            
------------------------------------------------------------------------
`min_c`   |    `-1.0`     | Minimum value for contour isovalue lines.   
          | (mode = ORB,  |                                             
          |   HFORB or    |                                             
          |    PROMPT)    |                                             
          |    `-0.1`     |                                             
          |  (MODE = QC,  |                                             
          |  INT or TOT)  |                                             
------------------------------------------------------------------------
`max_c`   |     `1.0`     | Minimum value for contour isovalue lines.   
          | (mode = ORB,  |                                             
          |   HFORB or    |                                             
          |    PROMPT)    |                                             
          |    `-0.1`     |                                             
          |  (MODE = QC,  |                                             
          |  INT or TOT)  |                                             
------------------------------------------------------------------------
`auto_c`  |     `no`      | If set `yes`, it defines `min_c` and `max_c`   
          |               | automatically according to the range of     
          |               | grid values.                                   
------------------------------------------------------------------------
`p_unit`  |     `angs`    | Length unity in the plot axes. If not set   
          |               | `angs` (Å), plot will be shown in atomic    
          |               | units (Bohr radius).                         
------------------------------------------------------------------------
`draw_atom`|    `plane`    | Which atoms will be shown on the plot.
           |               |                                             
           |               | Possible values:                            
           |               |                                             
           |               |`all`: All atoms are projected on the       
           |               |       plotted plane and shown on the plot. 
           |               |`plane`: Only atoms on the plotted plane and  
           |               |         0.1 Å or less away from the plotted  
           |               |         plane are shown on the plot.         
           |               |`none`: No atom is shown on the plot.        
 -----------------------------------------------------------------------
`draw_name`|     `yes`     | Whether atoms are represented on the plot   
           |               | by name (the labels given in $DATA section  
           |               | of GAMESS). If "no", atoms are represented  
           |               | by a black dot.                             
------------------------------------------------------------------------


#### 5.1.4 External File Settings

VARIABLE  | DEFAULT VALUE |  DESCRIPTION                                
------------------------------------------------------------------------
`save_png`|   -   (300)   | File name for saving the plot as a `*.png` 
          |               | image. The variable is optional and is not
          |               | used if mode = `PROMPT`. The second argument
          |               | must be separated by a blankspace and 
          |               | corresponds to the figure resolution in DPI
          |               | (dots per inch). The default value is 300.
------------------------------------------------------------------------
`save_eps`|   -   (300)   | File name for saving the plot as a `*.eps`    
          |               | image. The variable is optional and is not
          |               | used if mode = `PROMPT`. The second argument
          |               | must be separated by a blankspace and 
          |               | corresponds to the figure resolution in DPI
          |               | (dots per inch). The default value is 300.
------------------------------------------------------------------------
`save_txt`|       -       | File name for saving the plot in text format
          |               | (`*.txt`). The variable is optional and is not
          |               | used if mode = `PROMPT`. For information about
          |               | formatting, please see section 6.4.         
------------------------------------------------------------------------

Please check the file [TUTORIAL.md](TUTORIAL.md) for further information.

## 6. RESOURCES
------------

### 6.1 AVAILABLE COLORMAPS

Matplotlib provides a huge number of colormaps. They can be checked out [here](https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html).

The recommended colormaps for using in GPF-PLOT are shown in the image `examples/useful_colormaps.png`. My personal favorite is `coolwarm` for plotting orbitals and `seismic_r` for plotting densities. These are the default colormaps mentioned above.

*Hint*: all colormaps can be reversed by adding `_r` at the end of its name.


### 6.2 The gpfplot_figure SCRIPT

This script is an extra resource, used for generating high-quality figures intended for publishing in papers or books. The figure may contain more than one orbital or density plot.

Check [TUTORIAL.md](TUTORIAL.md) for information about how to use it.


## 7. FINAL REMARKS
----------------

There are some improvements wanted for a newer version, but I do not know when I will be able to work on it. Although, any suggestions are welcome. Since the program has not been fully tested, it probably may have bugs.

You can send questions and suggestions to [Dr. David Sousa](mailto:david.sousarj@yahoo.com.br). I will be glad to answer any questions and receive any suggestions on improving the program.
