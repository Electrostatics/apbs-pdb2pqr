SHELL=/bin/sh
EXE=${AMBERHOME}/exe/mm_pbsa.pl

#test.all: test.GenerateSnapshots test.MMPBSA_Stability test.MMPBSA_Binding \
#          test.MMPBSA_Decomp_Residue test.MMPBSA_Decomp_Pair

test.all: test.GenerateSnapshots test.MMPBSA_Stability test.MMPBSA_Binding

test.GenerateSnapshots:
	cd 01_GenerateSnapshots; $(EXE) mm_pbsa.in > mm_pbsa.log 2>&1

test.MMPBSA_Stability:
	cd 02_MMPBSA_Stability; $(EXE) mm_pbsa.in > mm_pbsa.log 2>&1

test.MMPBSA_Binding:
	cd 03_MMPBSA_Binding; $(EXE) mm_pbsa.in > mm_pbsa.log 2>&1

#test.MMPBSA_Nmode:
#	cd 04_MMPBSA_Nmode/; $(EXE) mm_pbsa.in > mm_pbsa.log 2>&1

#test.MMPBSA_Decomp_Residue:
#	cd 05_MMPBSA_Decomp_Residue/; $(EXE) mm_pbsa.in > mm_pbsa.log 2>&1

#test.MMPBSA_Decomp_Pair:
#	cd 06_MMPBSA_Decomp_Pair; $(EXE) mm_pbsa.in > mm_pbsa.log 2>&1
