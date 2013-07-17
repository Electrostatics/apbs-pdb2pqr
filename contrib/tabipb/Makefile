# makefile for bimpb
FC=ifort
CFLAGS=-I. -L.

ifndef DEBUG
endif
	
ifeq ($(shell uname -n), akira)
        FC=gfortran
	CFLAGS+=-O3
else
	CFLAGS+=-fast
endif

OBJS = stub.o var_modules.o treecode3d_pb.o readin.o dgmres_dep.o
NAME = tabipb
EXECUTABLE = $(NAME)
LIBNAME = lib$(EXECUTABLE)

.PHONY: all
all: tabipb

$(LIBNAME).a: $(OBJS)
	ar cr $@ $?

$(EXECUTABLE): $(LIBNAME).a $(EXECUTABLE).o
	$(FC) -o$@ $(EXECUTABLE).o -l$(NAME) $(CFLAGS)

%.o:%.f*
	$(FC) $(CFLAGS) -c $<

%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean: 
	-rm *.o *.mod $(EXECUTABLE) $(LIBNAME).* 
