# makefile for bimpb
F90= ifort
flag= -fast -c  

ifeq ($(shell uname -n), akira)
        F90=gfortran
	flag=-O3 -c
#        FCFLAGS+=-std=legacy -ffree-form -ffree-line-length-200
endif

tabipb.out: var_modules.o treecode3d_pb.o main.o readin.o dgmres_dep.o
	$(F90)  -o tabipb.out *.o 
var_modules.o:	var_modules.f90
	$(F90) $(flag) var_modules.f90
main.o:		main.f90
	$(F90) $(flag) main.f90
readin.o:	readin.f90
	$(F90) $(flag) readin.f90
treecode3d_pb.o:	treecode3d_pb.f
	$(F90) $(flag) treecode3d_pb.f
dgmres_dep.o:	dgmres_dep.f
	$(F90) $(flag) dgmres_dep.f
molecule.mod:   var_modules.f90
	$(F90) $(flag) var_modules.f90
comdata.mod:   var_modules.f90
	$(F90) $(flag) var_modules.f90
bicg.mod:   	var_modules.f90
	$(F90) $(flag) var_modules.f90
treecode.mod:	var_modules.f90
	$(F90) $(flag) var_modules.f90
treecode3d_procedures.mod:	treecode3d_pb.f
	$(F90) $(flag) treecode3d_pb.f
clean: 
	rm *.o *.mod tabipb.out 
