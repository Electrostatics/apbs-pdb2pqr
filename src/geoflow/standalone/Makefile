CFLAGS=-W -Wall -Wextra -I ../../../include/

ifdef DEBUG
	CFLAGS+=-ggdb3 -O0
else PROFILE
	CFLAGS+=-O1 -g -pg
else
	CFLAGS+=-O3
endif

OBJS=cpbconcz2.o modules.o pbsolvercz.o surfconcz.o solvationsub.o subfunction.o


.PHONY: all
all: osm.exe 

osm.exe: main.cpp $(OBJS)
	$(CXX) $(FCFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: ../%.cpp
	$(CXX) $(CFLAGS) -c $<

.PHONY:doc
doc:
	doxygen

.PHONY: clean
clean:
	-rm -rf *.o *.mod *.exe doc

