#! /usr/bin/env python3.8
# 
# GPFPlot Figure Generation
# Created by David W. O. de Sousa, david.sousarj@yahoo.com.br
# Version 0.2.1, February 2025.
#
# Imports and functions ################################################
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

# Uncomment the lines below if gpfplot_figure.py is not in the GPFPlot 
# parent directory. Edit the variable GPFPATH adequately.
#GPFPATH = "/home/david/bin/gpfplot-0.2.1"
#import sys
#sys.path.append(GPFPATH)

from core.parse_input import *
from core.parse_vb import parse_geom
from core.plot_utils import Plot_Function, Plot_Settings
from core.operations import RAng

# Read variables #######################################################
figure_name = "C2H4_sigmaCC"     

n_cols = 3

n_rows = 1

# List of *.txt and *.gpfplot files to open
p_list = ["example_fig1", "example_fig2", "example_fig3"]

# Title of every plot
titles = ["\\rho^{INT}", "\\rho^{QC}", "\\rho^{TOT}"]

# Superior Title 
suptitle = "C_2 H_4\ \sigma(C-C)\ Density\ Partitioning"

# print size of each graph in inches
p_size = 3 

# Dots per inch (resolution)
DPI = 300

# Generate EPS vector
EPS = True

# Create Figure ########################################################
fig = plt.figure(figsize=(1 + p_size*n_cols, p_size*n_rows))
gs = gridspec.GridSpec(n_rows, n_cols)

axes = []
for i in range(n_rows):
	for j in range(n_cols):
		axes.append( plt.subplot(gs[i,j]) )

suptitle = r"$\mathdefault{" + suptitle + "}$"
titles=[ r"$\mathdefault{" + i + "}$" for i in titles]

# For each p_list object ###############################################
for i in range(len(p_list)):
	# Parse text grids 
	f = open(p_list[i]+".txt")
	ln = f.readline()
	ln = f.readline()
	out_file = ln[23:-1]
	ln = f.readline()
	#dens_file = ln[23:-1]
	ln = f.readline()
	gpf_file = ln[23:-1]
	ln = f.readline()
	mode = ln[2:-1]
	ln = f.readline()
	plane = ln[9:11]
	ln = f.readline()
	Xrng = ln[11:-1]
	ln = f.readline()
	Yrng = ln[11:-1]
	ln = f.readline()
	offset = float( ln[11:-1] )
	ln = f.readline()
	gridp = int( ln[13:].split("x")[0] )
	f.close()

	Xlim = tuple(map(float, Xrng.replace("(","").replace(")","").split(",")))
	Ylim = tuple(map(float, Yrng.replace("(","").replace(")","").split(",")))
	XX = np.linspace(Xlim[0], Xlim[1], gridp)
	YY = np.linspace(Ylim[0], Ylim[1], gridp)
	X, Y = np.meshgrid(XX, YY)

	PSI = np.loadtxt(p_list[i]+".txt", skiprows=10)

	mode0 = mode.split(" ")[0]
	gpf_text = open(gpf_file,'r').read().split('\n')
	c_fill, color_f, c_lines, color_l, c_label,\
        min_c, max_c, nconts, auto_c,\
	    p_unit, draw_atom, draw_name = parse3(gpf_text, mode0)
	
	out_text = open(out_file,'r').readlines()
	NATOMS, CHARGE, ATOMS, R = parse_geom(out_text)

	# in Plot_Settings grid enters in bohrs
	X /= RAng; Y /= RAng

	Atom_posx, Atom_posy, Atom_labl,\
    label_x, label_y, X, Y = Plot_Settings(ATOMS, NATOMS, R,
                                           draw_atom, plane, offset,
                                           p_unit, X, Y)

	Plot_Function(axes[i], c_fill, color_f, c_lines, color_l, c_label,
                  min_c, max_c, nconts, auto_c, p_unit,
                  draw_atom, draw_name, Atom_posx, Atom_posy,
                  Atom_labl, label_x, label_y, X, Y, PSI, 1)

	axes[i].text(np.median(XX), Ylim[1]-0.1, titles[i],
                 verticalalignment='top',
                 horizontalalignment='center',
                 fontsize=14)

	# for all rows except the last one
	if i < (n_cols)*(n_rows-1):
		plt.setp(axes[i].get_xticklabels(), visible=False)
		axes[i].set_xlabel("")

	# for all columns except the first one
	if i not in list(range(0, n_cols*n_rows, n_cols)):
		plt.setp(axes[i].get_yticklabels(), visible=False)
		axes[i].set_ylabel("")

	print( p_list[i]+" PLOTTED." )
#end
fig.suptitle(suptitle, fontsize=16)
plt.tight_layout()

# Manual adjustments ###################################################

fig.subplots_adjust(left=0.08)
fig.subplots_adjust(bottom=0.08)
fig.subplots_adjust(right=0.96)
fig.subplots_adjust(top=0.92)
fig.subplots_adjust(wspace=0.16)
fig.subplots_adjust(hspace=0.17)

########################################################################

fig.savefig("%s.png" %figure_name, dpi=300)
fig.savefig("%s.eps" %figure_name, dpi=300)
plt.show()


