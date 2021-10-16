#! /usr/bin/python3.8
#
# GPF-Plot
# Created by David W. O. de Sousa, david.sousarj@yahoo.com.br
# Version 0.2, October 2021.
__version__ = '0.2'
__author__ = "David W. O. de Sousa"
#
# NEEDED IMPROVEMENTS ##################################################
#TODO: fix symmetry problem in parsing basis set

##TODO: allow changing parse3 variables in prompt mode
##TODO: Read orbitals from VB2000 .guess files
##TODO: allow non-cartesian planes
##TODO: options KIN_QC, KIN_INT, KIN_TOT

# Imports and functions ################################################
import sys
import matplotlib.pyplot as plt
from core.parse_input import *
from core.parse_vb import *
from core.operations import *
from core.plot_utils import *

def gen_txt(out_file, dens_file, gpf_file, mode,
            plane, Xlim, Ylim, offset, gridp,
            txt, PSI):
	import os
	from numpy import savetxt
	from datetime import datetime

	offplane = 'xyz'.replace( plane[0], "").replace( plane[1], "")
	header ="""\
Generated from GPF-PLOT version {0} at {1:%Y-%m-%d %H:%M:%S}
GAMESS/VB2000 file:  {2}
Density Matrix file: {3}
GPF-PLOT input file: {4}
{5}
plane: {6}
{7} range:  {8}
{9} range:  {10}
{11} offset: {12}
gridpoints: {13}x{13}""".format(__version__, datetime.now(),
                                os.path.abspath(out_file),
                                os.path.abspath(dens_file),
                                os.path.abspath(gpf_file), mode,
                                plane, plane[0], Xlim, plane[1], Ylim,
                                offplane, offset, gridp)
	txt += ".txt"
	np.savetxt(txt, PSI, fmt="%11.6f", delimiter="", header=header) 
#end

# Splash screen ########################################################
print("""\
 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 -                                                                     -
 *               G    P    F    -    P    L    O    T                  *
 -                                                                     -
 *                            Version 0.2                              *
 -                            October 2021                             -
 *                                                                     *
 -                      Created by David Sousa                         -
 *                Inspired by G. N. Freitas' DENSPLOT                  *
 -             Also inspired by Brett Bode's wxMacMolPlt               -
 *                        Powered by Python                            *
 -                                                                     -
 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
""")

# INPUT READING ########################################################
try:
	arg1 = sys.argv[1]
except:
	print("Please specify input (*.gpfplot) file.")
	exit(1)
 
infile = open(arg1+".gpfplot", 'r').read().split('\n')

# input file and mode
out_file, dens_file, mode = parse1(infile)
mode0 = mode.split(" ")[0]

# Print input options
print()
if out_file == "":
	print(" ERROR: GAMESS / VB2000 input file not specified.")
	exit(1)
else:
	print(" GAMESS / VB2000 input file: " + out_file)

if dens_file == "":
	print(" GPF-EP Density matrix file: not specified.")
	if mode0 not in ("HFORB","ORB","PROMPT"):
		print(" ERROR: Density matrix file is required for mode "+mode0)
		exit(1)
	else:
		print(" Only ORB and HFORB options are available.")
else:
	print(" GPF-EP Density matrix file: {0}".format(dens_file))	

if mode0 not in MODES:
	print(" ERROR: Invalid mode.")
	exit(1)
else:
	print(" Mode: " + mode)
print()

# plane to plot
plane, gridp, Xlim, Ylim, offset = parse2(infile)

# graphic options input
c_fill, color_f, c_lines, color_l, c_label,\
min_c, max_c, nconts, auto_c,\
p_unit, draw_atom, draw_name = parse3(infile, mode0)

save_png, save_eps, dpipng, dpieps, save_txt = parse4(infile)
########################################################################
# PROGRAM START

# Parse input file
print(" Loading input files...")
jobfile = open(out_file,'r')
jobtext = jobfile.readlines()

NATOMS, CHARGE, ATOMS, R = parse_geom(jobtext)
L, M, N, A, IATOM, IATOM2, CC = parse_basis(jobtext)

if len(IATOM2) != NATOMS: #  symmetry issues with basis set
	L, M, N, A,\
    IATOM, CC = fix_basis(L, M, N, A, IATOM, IATOM2, ATOMS, CC)

NBF, GL, C, hf, g99 = parse_gpf(jobtext)
CHF = parse_hf(jobtext, NBF)

if dens_file != "":
	jobfile2 = open(dens_file,'r')
	jobtext2 = jobfile2.readlines()
	if g99: del(GL[-1])
	OVERLAPS = np.array( parse_overlaps(jobtext, hf, GL) )
	D_MATRIX = np.array( parse_dmatrix(jobtext2, GL) )
print()

# Calculate orbitals
print(" Calculating orbitals...")
X1,Y1,Pxyz,Ri2=calcgrid(plane,Xlim,Ylim,gridp,offset, R, L,M,N, IATOM)
NN=NNorm(L,M,N,A,CC,NBF)

AO_List = calc_AOs(NBF, A, CC, Ri2, NN, Pxyz, gridp)
NC, NHF = len(C), len(CHF)

if mode0 == "PROMPT":
	ORB_List= calc_all(NC, NBF, AO_List, C, gridp)
	HFO_List= calc_all(NHF, NBF, AO_List, CHF, gridp)
	print(" There are {0} HF and {1} VB orbitals available.".\
          format(NHF,NC))
else:
	modetxt=mode
	mode=mode.split()
	# mode ORB
	if mode[0] == 'ORB':
		norb = int(mode[1]) - 1	
		PSI = calc_orb(NBF, AO_List, C, norb, gridp)
	# mode HFORB
	elif mode[0] == 'HFORB':
		norb = int(mode[1]) - 1
		PSI = calc_orb(NBF, AO_List, CHF, norb, gridp)
	# mode QC, INT, TOT
	elif mode[0] in DENS:
		N = list(map(int, mode[1:]))

		# calculate just the needed orbitals
		ORB_List = []
		for a in range(NC):
			ORB_List.append([])
			if a+1 in N:
				ORB_List[a] = calc_orb(NBF, AO_List, C, a, gridp)

		if mode[0] == "QC":
			PSI = calc_QC(N, ORB_List, OVERLAPS, D_MATRIX, gridp)
		elif mode[0] == "INT":
			PSI = calc_INT(N, ORB_List, OVERLAPS, D_MATRIX, gridp)
		elif mode[0] == "TOT":
			PSI = calc_TOT(N, ORB_List, D_MATRIX, gridp)
	#end
#end
print()

# Graphic setup
print(" Setting up plot...")
Atom_posx, Atom_posy, Atom_labl,\
label_x, label_y, X1, Y1 = Plot_Settings(ATOMS, NATOMS, R, draw_atom,
                                         plane, offset, p_unit, X1, Y1)
print()


########################################################################
# Main Program
# Prompt mode
if mode == "PROMPT":
	while True:
		mode1 = input("gpfplot [{0}]> ".format(out_file[:-4]))	
		mode1 = mode1.split()

		if mode1 == []: # just pressed enter
			continue

		elif mode1[0] == "exit":
			exit()

		else:
			# Create Plot Window
			fig = plt.figure()

			# mode ORB
			if mode1[0] == 'ORB':
				norb = int( mode1[1] ) - 1
				fig.canvas.set_window_title("{0} ORB {1}".\
                                          format(out_file[:-4], norb+1))
				PSI = ORB_List[ norb ]

			# mode HFORB
			elif mode1[0] == 'HFORB':
				norb = int( mode1[1] ) - 1
				fig.canvas.set_window_title("{0} HFORB {1}".\
                                          format(out_file[:-4], norb+1))
				PSI = HFO_List[ norb ]

			# mode QC, INT, TOT
			elif mode1[0] in DENS:
				if dens_file == "":
					print("Error. No density matrix file defined.")
					plt.close()
					continue			

				N = list( map(int, mode1[1:]) )
				fig.canvas.set_window_title("{0} {2} {1}".\
                                            format(out_file[:-4],
                                            " ".join(str(N)), mode1[0]))
				if mode1[0] == "QC":
					PSI = calc_QC(N, ORB_List, OVERLAPS, D_MATRIX, gridp)
				elif mode1[0] == "INT":
					PSI = calc_INT(N, ORB_List, OVERLAPS, D_MATRIX, gridp)
				elif mode1[0] == "TOT":
					PSI = calc_TOT(N, ORB_List, D_MATRIX, gridp)
			else:
				print("Invalid mode! Try again.")
				plt.close()
				continue

			Plot_Function(plt, c_fill, color_f, c_lines, color_l,
                          c_label, min_c, max_c, nconts, auto_c, p_unit,
                          draw_atom, draw_name, Atom_posx, Atom_posy,
                          Atom_labl, label_x, label_y, X1, Y1, PSI, 0)
			plt.show()
			fig.clf()

# Single plot mode
else:
	fig = plt.figure()
	# mode ORB
	if mode[0] == 'ORB':
		fig.canvas.set_window_title("{0} ORB {1}".\
                                    format(out_file[:-4], norb+1))
	# mode HFORB
	elif mode[0] == 'HFORB':
		fig.canvas.set_window_title("{0} HFORB {1}".\
                                    format(out_file[:-4], norb+1))
	# mode QC, INT, TOT
	elif mode[0] in DENS:
		fig.canvas.set_window_title("{0} {2} {1}".\
                                    format(out_file[:-4],
                                           " ".join(str(N)), mode[0]))
	# error
	else:
		print("Invalid mode! Try again.")
		plt.close()
		exit(1)

	Plot_Function(plt, c_fill, color_f, c_lines, color_l, c_label,
                  min_c, max_c, nconts, auto_c, p_unit,
                  draw_atom, draw_name, Atom_posx, Atom_posy,
                  Atom_labl, label_x, label_y, X1, Y1, PSI, 0)
	print(" "+" ".join(mode)+" plotted.")

	if save_png != "":
		fig.savefig(save_png+".png", dpi=dpipng)
		print(" {0}.png plotted.".format(save_png))

	if save_eps != "":
		fig.savefig(save_eps+".eps", dpi=dpieps)
		print(" {0}.eps plotted.".format(save_eps))

	if save_txt != "":
		gen_txt(out_file, dens_file, arg1+".gpfplot", modetxt,
                plane, Xlim, Ylim, offset, gridp, save_txt, PSI)
		print(" {0}.txt saved.".format(save_txt))

	plt.show()
# The End!
