#! /usr/bin/env python3
# GPFPlot Library for Orbital and Density Calculations
# Last modified: 2021-09-30
#
# Imports and functions ################################################
import numpy as np

# Global variables and functions #######################################

DEBUG=False # Print text results

pi=np.pi
sqrt=np.sqrt
exp=np.exp
f=np.math.factorial
M=np.multiply
RAng=0.529177249 # angstrom / bohr
floor = lambda x: int(np.math.floor(x))

# double factorial
def df(n):
	if n <= 1:
		return 1
	elif n<=3:
		return n
	else:
		return n*df(n-2)

# Normalization for a uncontracted gaussian
def Norm(L,M,N,A):
	NN = (2*A/pi)**0.75
	if L!=0 and M!=0 and N!=0:
		NN *= 2**(L+M+N)
		NN *= sqrt(A**(L+M+N)/(df(L)*df(M)*df(N)))
	return NN

# Normalization for a contracted gaussian
def CNorm(L,M,N,A,CC,K):
	NN=0
	for i in range(K):
		for j in range(K):
			NN += CC[i]*CC[j]/(A[i]+A[j])**(L+M+N+1.5)
		#end
	#end
	NN *= df(2*L-1)*df(2*M-1)*df(2*N-1)*pi**0.75
	NN /= 2**(L+M+N)
	NN  = 1.0/sqrt(NN)
	return NN

# Normalization for each AO
def NNorm(L,M,N,A,CC,NBF):
	NN = np.zeros( NBF)
	for i in range(NBF):
		NN[i]=CNorm(L[i],M[i],N[i],A[i],CC[i],len(CC[i]))
	#end
	return NN

# Generate grid for each atom: Ri^2 and polinomial part of orbital
def calcgrid(plane,Xlim,Ylim,gridp,offset,R, L,M,N, IATOM):
	# for plotting
	Xrange = np.linspace(*Xlim,num=gridp)
	Yrange = np.linspace(*Ylim,num=gridp)
	Xrange = Xrange/RAng # convert to bohr!
	Yrange = Yrange/RAng # convert to bohr!
	X1, Y1 = np.meshgrid(Xrange,Yrange)

	# Construct grids
	if plane == "xy":
		X=X1; Y=Y1; Z=offset
	elif plane == "yx":
		X=Y1; Y=X1; Z=offset
	elif plane == "yz":
		X=offset; Y=X1; Z=Y1
	elif plane == "zy":
		X=offset; Y=Y1; Z=X1
	elif plane == "xz":
		X=X1; Y=offset; Z=Y1
	elif plane == "zx":
		X=Y1; Y=offset; Z=X1
	else:
		print(" Error: Plane not valid.")
		exit(1)
	#endif

	# number of atoms
	NBF  = len(IATOM)

	Ri2 = np.zeros([NBF,gridp,gridp])
	P   = np.zeros([NBF,gridp,gridp])

	for i in range(NBF):
		P[i]   = (X - R[0][IATOM[i]])**L[i] * (Y - R[1][IATOM[i]])**M[i] * (Z - R[2][IATOM[i]])**N[i] 
		Ri2[i] = (X - R[0][IATOM[i]])**2 + (Y - R[1][IATOM[i]])**2 + (Z - R[2][IATOM[i]])**2

	# DEBUG ############################################################
	if DEBUG:
		P.tofile("tempfiles/Pxyz.txt",sep=" ", format="%10.6f")
		Ri2.tofile("tempfiles/Ri2.txt",sep=" ", format="%10.6f")
	####################################################################
	return X1,Y1, P, Ri2
		
# Calculate atomic orbitals
def calc_AOs(NBF, A, CC, Ri2, NN, Pxyz, gridp):
	AO_List=[]

	for i in range(NBF):
		sumj=np.zeros( [ gridp,gridp ] )
		for j in range(len(CC[i])):
			sumj=sumj+CC[i][j]*exp(-A[i][j]*Ri2[i]) 
		#end
		AO=NN[i]*Pxyz[i]*sumj #
	# DEBUG ############################################################
		if DEBUG:
			AO.tofile("tempfiles/AO%d.txt" %i,sep=" ", format="%10.6f")
	####################################################################
		AO_List.append( AO )
	#end
	return AO_List

# Calculate orbital in terms of AOs
def calc_orb(NBF, AO_List, C, norb, gridp):
	# define wavefunction grid
	PSI = np.zeros( [ gridp,gridp ] )
	for i in range(NBF):
		PSI=PSI+C[norb][i]*AO_List[i]
	#end
	return PSI

# Calculate all orbitals
def calc_all(N, NBF, AO_List, C, gridp):
	ORB_List=[]
	for norb in range(N):
		PSI= calc_orb(NBF, AO_List, C, norb, gridp)
		# DEBUG ########################################################
		if DEBUG:
			PSI.tofile("tempfiles/PSI%d.txt" %norb,sep=" ", format="%10.6f")
		################################################################
		ORB_List.append( PSI )
	return ORB_List

# Calculate total density of a set of orbitals
def calc_TOT(N, PHI, P, gridp):
	# N   = Set of orbitals (by index, starting at 1) 
	# PHI = list of orbital grids
	# P   = one-electron density matrix in the orbital basis
	RHO = np.zeros( [gridp, gridp] )
	for i in N:
		r = i -1
		for j in N:
			s = j - 1
			RHO = RHO + M(PHI[r], PHI[s])*P[r][s]			
	return RHO

# Calculate quasi-classical density of a set of orbitals
def calc_QC(N, PHI, S, P, gridp):
	# N   = Set of orbitals (by index, starting at 1) 
	# PHI = list of orbital grids
	# S   = overlap matrix
	# P   = one-electron density matrix in the orbital basis
	RHO = np.zeros( [gridp, gridp] )
	for i in N:
		r = i - 1
		for j in N:
			s = j - 1
			RHO = RHO + M(PHI[r], PHI[r])*P[r][s]*S[r][s]			
	return RHO

# Calculate interference density of a set of orbitals
def calc_INT(N, PHI, S, P, gridp):
	# N   = Set of orbitals (by index, starting at 1) 
	# PHI = list of orbital grids
	# S   = overlap matrix
	# P   = one-electron density matrix in the orbital basis
	RHO = np.zeros( [gridp, gridp] )
	for i in N:
		r = i - 1
		PHIr2 = M(PHI[r], PHI[r])
		for j in N:
			s = j - 1
			if r != s:
				PHIs2    = M(PHI[s], PHI[s])
				PHIrPHIs = M(PHI[r], PHI[s])
				RHO = RHO + (PHIrPHIs -.5*S[r][s]*(PHIr2 + PHIs2))*P[r][s]			
	return RHO
