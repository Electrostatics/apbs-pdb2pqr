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

OBJS = var_modules.o treecode3d_pb.o readin.o dgmres_dep.o
NAME = tabipb
EXECUTABLE = $(NAME)
LIBNAME = lib$(EXECUTABLE)

.PHONY: all
all: $(LIBNAME).a

$(LIBNAME).a: $(OBJS)
	ar cr $@ $?

#$(EXECUTABLE): $(LIBNAME).a $(EXECUTABLE).o
#	$(FC) -o$@ $(EXECUTABLE).o -l$(NAME) $(CFLAGS)

%.o:%.f*
	$(FC) -c $< $(CFLAGS)

%.o:%.c
	$(CC) -c $< $(CFLAGS)

clean: 
	-rm *.o *.mod $(EXECUTABLE) $(LIBNAME).* 
