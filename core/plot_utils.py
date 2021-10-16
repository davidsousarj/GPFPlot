#! /usr/bin/env python3
# GPFPlot Library for Plotting Utilities
# Last modified: 2020-07-09
#
import numpy as np
from matplotlib.pyplot import colorbar
import matplotlib.cm as cm
from core.operations import RAng

def Plot_Settings(ATOMS, NATOMS, R,
                  draw_atom, plane, offset,
                  p_unit, X1, Y1):
	AxisIndex = {'x':0, 'y':1, 'z':2}
	if draw_atom == 'all':
		Atom_posx = np.array( R[ AxisIndex[ plane[0] ] ][:] )
		Atom_posy = np.array( R[ AxisIndex[ plane[1] ] ][:] )
		Atom_labl = ATOMS
	elif draw_atom == 'plane':
		offplane = 'xyz'.replace( plane[0], "").replace( plane[1], "")
		Atom_posx = []
		Atom_posy = []
		Atom_labl = []
		for i in range(NATOMS):
			if np.abs( R[ AxisIndex[ offplane ] ][i] - offset ) <= 0.2:
				Atom_posx.append( R[ AxisIndex[ plane[0] ] ][i] )
				Atom_posy.append( R[ AxisIndex[ plane[1] ] ][i] )
				Atom_labl.append( ATOMS[i] )
		Atom_posx = np.array( Atom_posx )
		Atom_posy = np.array( Atom_posy )
	elif draw_atom == 'none':
		Atom_posx = []
		Atom_posy = []
		Atom_labl = []
		
	if p_unit == 'angs':
		X1 *= RAng
		Y1 *= RAng
		Atom_posx *= RAng
		Atom_posy *= RAng

	label_0=' / Å'
	label_x =plane[0]+label_0
	label_y =plane[1]+label_0

	return Atom_posx, Atom_posy, Atom_labl, label_x, label_y, X1, Y1

def Plot_Function(plt, c_fill, color_f, c_lines, color_l, c_label,
                  min_c, max_c, nconts, auto_c, p_unit,
                  draw_atom, draw_name, Atom_posx, Atom_posy,
                  Atom_labl, label_x, label_y, X1, Y1, PSI, gs):
	TextProp = {'horizontalalignment':'center',
                'verticalalignment':'center'}
	LabelProp= {'size':18}

	if not auto_c:
		levls=np.linspace(min_c, max_c, nconts)
		cticks=np.linspace(min_c, max_c, 11)
	else:
		lim  = np.abs( np.amax(PSI) )
		lim2 = np.abs( np.amin(PSI) )
		lim = max(lim, lim2)
		if lim == 0.: lim = 1.
		levls=np.linspace(-lim, lim, nconts)
		cticks=np.linspace(-lim, lim, 11)
	if c_fill:
		cf=plt.contourf(X1,Y1,PSI, cmap=eval(color_f), levels=levls)
		# Colorbar syntax is different between GridSpec and simple plot.
		if gs == 0:
			plt.colorbar(ticks=cticks, pad=0)
		else:
			colorbar(cf, ax=plt, ticks=cticks, pad=0)
	if c_lines:
		if c_fill:
			width = 0.2
		else:
			width = 0.5
		if color_l == "":
			CS = plt.contour(X1,Y1,PSI,
                             colors='k',
                             levels=levls,
                             linewidths=width) 
		else:
			CS = plt.contour(X1,Y1,PSI,
                             linewidths=width,
                             cmap=color_l,
                             levels=levls)
		if c_label:
			plt.clabel(CS, inline=1, fontsize=10)
	if gs == 0:
		plt.xlabel(label_x)
		plt.ylabel(label_y)
	else:
		plt.set_xlabel(label_x)
		plt.set_ylabel(label_y)

	# plot atom positions
	if draw_atom != 'none':
		for i in range(len(Atom_posx)):
			if draw_name:
				plt.text(Atom_posx[i], Atom_posy[i],
                         Atom_labl[i], **TextProp, **LabelProp)
			else:
				plt.plot(Atom_posx[i], Atom_posy[i], 'ko')
#end
