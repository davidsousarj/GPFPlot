#! /usr/bin/env python3
# GPFPlot Library for Parsing Input Files
# Last modified: 2020-07-09

from math import ceil
from core.parse_input import getLine
from core.operations import RAng
# Parameters ###########################################################
# basis functions per shell type
ShellSize = {"S" : 1, "P" : 3, "L" : 4, "D" : 6, "F" : 10, "G" : 15,
             "H": 21, "I" : 28 }

# L, M,N vector templates according with shell type
LMN = [\
[[0], [0], [0]],
[[1,0,0], [0,1,0], [0,0,1]],
[[0,1,0,0], [0,0,1,0], [0,0,0,1]],
[[2,0,0,1,1,0], [0,2,0,1,0,1], [0,0,2,0,1,1]],
[[3,0,0,2,2,1,0,1,0,1], [0,3,0,1,0,2,2,0,1,1], [0,0,3,0,1,0,1,2,2,1]],
[[4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1],
 [0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1],
 [0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2]],
[[5, 0, 0, 4, 4, 1, 0, 1, 0, 3, 3, 2, 0, 2, 0, 3, 1, 1, 2, 2, 1],
 [0, 5, 0, 1, 0, 4, 4, 0, 1, 2, 0, 3, 3, 0, 2, 1, 3, 1, 2, 1, 2],
 [0, 0, 5, 0, 1, 0, 1, 4, 4, 0, 2, 0, 2, 3, 3, 1, 1, 3, 1, 2, 2]],
[[6, 0, 0, 5, 5, 1, 0, 1, 0, 4, 4, 2, 0, 2,
  0, 4, 1, 1, 3, 3, 0, 3, 3, 2, 1, 2, 1, 2],
 [0, 6, 0, 1, 0, 5, 5, 0, 1, 2, 0, 4, 4, 0,
  2, 1, 4, 1, 3, 0, 3, 2, 1, 3, 3, 1, 2, 2],
 [0, 0, 6, 0, 1, 0, 1, 5, 5, 0, 2, 0, 2, 4,
  4, 1, 1, 4, 0, 3, 3, 1, 2, 1, 2, 3, 3, 2]]]

# index of shell types in LMN array
ShellIndex = {"S" : 0, "P" : 1, "L" : 2, "D" : 3,
              "F" : 4, "G" : 5, "H" : 6, "I" : 7}

# Functions ############################################################
def parse_geom(jobtext):
	""" Get geometric parameters.""" 
	line = getLine(jobtext," TOTAL NUMBER OF ATOMS                    ")
	NATOMS = int( jobtext[line].split()[-1] )

	line = getLine(jobtext," ATOM      ATOMIC")
	unit = jobtext[line].split()[-1][1:-1]

	LATOMS = []
	CHARGE = []
	COORD  = [[],[],[]]
	for i in range(NATOMS):
		LATOMS.append( jobtext[line + 2 + i].split()[0] )
		CHARGE.append( float( jobtext[line + 2 + i].split()[1] ))
		COORD[0].append( float( jobtext[line + 2 + i].split()[2] ))
		COORD[1].append( float( jobtext[line + 2 + i].split()[3] ))
		COORD[2].append( float( jobtext[line + 2 + i].split()[4] ))
	# converting angs into bohrs (if applicable)
	if unit == 'ANGS':
		for i in range(3):
			for j in range(NATOMS):
				COORD[i][j] /= RAng

	return NATOMS, CHARGE, LATOMS, COORD

def get_symm(jobtext):
	"""
	Get point group symmetry of the molecule
	in order to fix the basis set."""
	ln = getLine(jobtext, "THE POINT GROUP OF THE MOLECULE IS")
	pg = jobtext[ln].split()[-1]
	return pg

def parse_basis(jobtext):
	""" Get basis set parameters. """
	L = []; M = []; N = [] 	# exponents for x, y, z
	A = []             		# exponent of (r-R)^2
	IATOM = []          	# which atom (index)
	CC = []             	# contraction coefficients
	IATOM2 = []             # which atom (name) - fix symmetry problem

	li = getLine(jobtext, "ATOMIC BASIS SET")+7
	lf = getLine(jobtext, "TOTAL NUMBER OF BASIS")
	atomcount = 0
	currentshell = 0

	for ln in range(li,lf):
		t = jobtext[ln].split()

		# pass blank line
		if not len(t):
			continue

		# if contains any text, change current atom
		elif t[0].upper().isupper():
			atomcount += 1
			IATOM2.append(t[0])

		# line with numeric info
		elif t[0].isdigit():
			i  = ShellIndex[t[1]] # index of LMN exponents
			sz =  ShellSize[t[1]] # size of this shell

			# if it is a new shell
			if t[0] != currentshell:
				currentshell = t[0]
				# add corresponding values to vectors L,M,N
				L.extend(LMN[i][0])
				M.extend(LMN[i][1])
				N.extend(LMN[i][2])
				# add atom geometry "sz" times to R array
				IATOM.extend([atomcount-1]*sz)
				# add corresponding values to A,CC arrays "sz" times
				# python syntax [[x]]*y creates y virtual copies of [x]
				# if one [x] is changed, all of them will be as well
				A.extend([[float(t[3])]]*sz)
				CC.extend([[float(t[4])]]*sz)

			# continue last currentshell
			# add only to X[-1] but the "sz" last positions will be
			# changed, due to the python syntax [[x]]*y
			elif t[0] == currentshell:
				A [-1].append(float(t[3]))
				CC[-1].append(float(t[4]))

	# The End!Spit it all out!
	return L, M, N, A, IATOM, IATOM2, CC

def fix_basis(L, M, N, A, IATOM, LABELS, NEW_LABELS, CC):
	"""
	Fix basis set shown only for symmetry-
	unique atoms."""
	# LABELS is the current layout of the basis set vectors
	# it must be transformed to the layout given by NEW_LABELS

	# Calculate basis size of each atom
	atom_bs = []
	compare = -1
	for i in range(len(IATOM)):
		if IATOM[i] != compare:
			compare += 1
			atom_bs.append(i)
	atom_bs.append(len(IATOM))

	# Generate new basis set
	new_L = []; new_M = []; new_N = [];
	new_A = []; new_CC = [];
	for atom in NEW_LABELS:
		index = LABELS.index(atom)
		new_L.extend(L[atom_bs[index] : atom_bs[index+1]])		
		new_M.extend(M[atom_bs[index] : atom_bs[index+1]])
		new_N.extend(N[atom_bs[index] : atom_bs[index+1]])
		new_A.extend(A[atom_bs[index] : atom_bs[index+1]])
		new_CC.extend(CC[atom_bs[index] : atom_bs[index+1]])

	# Generate new IATOM
	new_IATOM = []
	for i in range(len(NEW_LABELS)):
		index = LABELS.index(NEW_LABELS[i])
		new_IATOM.extend([i]*(atom_bs[index+1] - atom_bs[index]))

	return new_L, new_M, new_N, new_A, new_IATOM, new_CC


def	parse_gpf(jobtext):
	line = getLine(jobtext, "NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCT")
	NBF = int(jobtext[line].split()[-1])
	line = getLine(jobtext, " NUMBER OF ELECTRONS                    ")
	NELS = int(jobtext[line].split()[-1])
	line = getLine(jobtext," Number of electron groups       =")
	NGROUPS = int(jobtext[line].split()[-1])
	line = getLine(jobtext," Num. of electron")
	EL_GROUP = list(map(int, jobtext[line].split()[3:]))
	line = getLine(jobtext," Num. of orbital")
	ORB_GROUP = list( map(int, jobtext[line].split()[3:]))
	NORBS = sum(ORB_GROUP)

	GROUP_LAYOUT = []
	count = 0
	for i in range(NGROUPS):
		GROUP_LAYOUT.append(list(range(count, ORB_GROUP[i] + count )))
		count += ORB_GROUP[i]

	line = getLine(jobtext," Method# ")
	hf_group = int(jobtext[line].split()[1])

	# If SPHER is set an extra group with Method#=99 is created and
	# the overlaps are not printed
	line = getLine(jobtext," Method# ")
	g99 = jobtext[line].split()[-1] == "99"

	lo = getLine(jobtext, "ORBITALS OF EACH ELECTRON")+6
	ncols=6
	space=6

	# Ns is the number of horizontal orbital sets found in jobtext
	# Nl is the arrangement of orbitals found in jobtext
	Ns = int(ceil(NORBS / float(ncols)))
	Nl = [ncols]*Ns
	if NORBS % ncols != 0: Nl[-1] = NORBS % ncols

	# Extract the orbital coefficients
	C   = []
	for i in range(Ns):
		for j in range(Nl[i]):
			C.append([])
			for k in range(NBF):
				orbs = jobtext[lo+k][16:].split()
				C[j + ncols*i].append(float(orbs[j])) 
		lo += NBF + space # pass to next horizontal set of orbitals
	return NBF, GROUP_LAYOUT, C, hf_group, g99

def	parse_hf(jobtext, NBF):
	try:
		ln = getLine(jobtext, "KEPT IN THE VARIATION SPACE IS",
                     none=False)
	except:
		ln = getLine(jobtext, "NUMBER OF CARTESIAN", none=False)
	NORBS = int( jobtext[ln].split()[-1] ) # number of orbitals
	try:
		lo = getLine(jobtext, "MOLECULAR ORBITALS")+6
	except:
		lo = getLine(jobtext, "EIGENVECTORS")+6

	ncols=5
	space=4
	# Ns is the number of horizontal orbital sets found in jobtext
	# Nl is the arrangement of orbitals found in jobtext
	Ns = int(ceil(NORBS / float(ncols)))
	Nl = [ncols]*Ns
	if NORBS % ncols != 0: Nl[-1] = NORBS % ncols
	# Extract the orbital coefficients
	C   = []
	for i in range(Ns):
		for j in range(Nl[i]):
			C.append([])
			for k in range(NBF):
				orbs = jobtext[lo+k][16:].split()
				C[j + ncols*i].append(float(orbs[j])) 
		lo += NBF + space # pass to next horizontal set of orbitals
	return C

def parse_overlaps(jobtext, hf, GL):
	NO=GL[-1][-1] + 1	# number of orbitals
	NG = len(GL)  		# number of groups

#	#OVERLAPS=[] #np.identity(NO)
	OVERLAPS=[]
	for i in range(NO):
		OVERLAPS.append([])
		for j in range(NO):
			if i == j:
				OVERLAPS[i].append(1.0)
			else:
				OVERLAPS[i].append(0.0)	

	if hf == 1: # no HF group
		start = 1
	else:
		start = 0
	for i in range(start, NG):
		line = getLine(jobtext,
                     "OVERLAP MATRIX OF VB ORBITALS FOR GROUP  {0:2d}".\
                     format(i+1)) + 5
		for j in range(len(GL[i]) - 1):
			for k in range(j+1):
				OVERLAPS[GL[i][j + 1]][GL[i][k]]\
                         = float(jobtext[line][4 + 7*k : 11 + 7*k])
				OVERLAPS[GL[i][k]][GL[i][j + 1]]\
                         = float(jobtext[line][4 + 7*k : 11 + 7*k])
			#end
			line+=1
		#end
	#end
	return OVERLAPS	

def parse_dmatrix(jobtext, GL):
	NO=GL[-1][-1] + 1	# number of orbitals
	NG = len(GL)  		# number of groups

	#DMATRIX= np.zeros([NO,NO])
	DMATRIX=[]
	for i in range(NO):
		DMATRIX.append([])
		for j in range(NO):
			DMATRIX[i].append(0.0)	

	for i in range(NG):
		line = getLine(jobtext,
                       "ONE-ELECTRON DENSITY MATRIX - GROUP   {:3d}".\
                       format(i+1)) + 3
		if len(GL[i]) <= 10:
			for j in range(len(GL[i])):
				for k in range(len(GL[i])):
					DMATRIX[GL[i][j]][GL[i][k]]\
                            = float(jobtext[line][3 + 13*k : 16 + 13*k])
				#end
				line += 1
			#end
		# different format for groups > 10 orbitals
		else:
			worklist = [10]*(len(GL[i])//10)
			worklist.append(len(GL[i]) % 10)
			for w in range(len(worklist)):
				for j in range(len(GL[i]) - 10*w):
					for k in range(worklist[w]):
						DMATRIX[GL[i][j + 10*w]][GL[i][k + 10*w]]\
		                    = float(jobtext[line][3 + 13*k : 16 + 13*k])
					#end
					line += 1
				#end
				line += 2
			#end
		#end			
	#end
	return DMATRIX
