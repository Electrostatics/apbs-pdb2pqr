FC=gfortran
CFLAGS=-W -Wall

ifdef DEBUG
	CFLAGS+=-ggdb3 -O0
else
	CFLAGS+=-O3
endif

FCFLAGS+=-std=legacy -ffree-form -ffree-line-length-200
LDFLAGS=-L. -lslatec -lc -lstdc++

OBJS=adi.f90 fish.f90 pbsolvercz.f90 solvationsub.f90 surfconcz.f90 fftpack.f90 hw3crt.f90 pbconcz2.f90 pois3d.f90 subfunction.f90

.PHONY: all
all: osm.exe 

osm.exe:   $(OBJS)
	$(FC) $(FCFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) $(CFLAGS) -c $<

%.o: %.c
	$(CXX) $(CFLAGS) -c $<


.PHONY: clean
clean:
	-rm -f *.o *.mod *.exe
