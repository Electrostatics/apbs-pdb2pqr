## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .cc .c .f .o
INCLUDES           = @inc_list@
LDFLAGS            = @LDFLAGS@ @FLIBS@
CCLD = $(F77)

noinst_PROGRAMS    = mgmesh dxmath
mgmesh_SOURCES     = mgmesh.c 
mgmesh_LDADD       = ../../src/aaa_lib/libapbs.la @lib_list@
dxmath_SOURCES     = dxmath.c 
dxmath_LDADD       = ../../src/aaa_lib/libapbs.la @lib_list@

