## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .cc .c .f .o
CCLD = $(F77) $(FFLAGS)

build_libraries = @apbs_build_libraries@

bindir = ${prefix}/bin
bin_PROGRAMS     = apbs
apbs_SOURCES     = main.c 
apbs_LDADD       = ${build_libraries}

# This is a hack I suspect, but needed to force apbs to be
# rebuilt anytime the make command is issued
# Problem before is that there is no dependency between source file
# recompilation and relinking the binary. There is probably
# a cleaner way to do this, I can't find it, however.
BUILT_SOURCES	 = rm.apbs

rm.apbs: Makefile
	-rm -f apbs

libdir = ${prefix}/lib
lib_LTLIBRARIES = libapbsmainroutines.la
libapbsmainroutines_la_SOURCES = routines.c
libapbsmainroutines_la_LDFLAGS = -version-info ${APBS_VERSION}

INCLUDES = @inc_list@
AM_CFLAGS= @wall_pedantic@
AM_LDFLAGS       = @LDFLAGS@ @FFLAGS@
 
includedir = ${prefix}/include/apbs
include_HEADERS = routines.h
	
install-data-local:
	-cp -p ${top_apbsdir}/bin/ApbsClient.py ${prefix}/bin
	-chmod 755 ${prefix}/bin/ApbsClient.py
	-rm -f ${prefix}/bin/wsdl2py 