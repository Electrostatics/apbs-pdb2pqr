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
AM_CFLAGS = @CFLAGS@
AM_LDFLAGS = @LDFLAGS@ 

noinst_PROGRAMS          = mgmesh dxmath uhbd_asc2bin smooth value dx2mol benchmark similarity analysis
mgmesh_SOURCES           = mgmesh.c 
mgmesh_LDADD             = ../../src/aaa_lib/libapbs.la @lib_list@
mgmesh_LDFLAGS           = ../../src/aaa_lib/libapbs.la @lib_list@
mgmesh_CCLD              = $(CC)
smooth_SOURCES           = smooth.c 
smooth_LDADD             = ../../src/aaa_lib/libapbs.la @lib_list@
smooth_LDFLAGS           = ../../src/aaa_lib/libapbs.la @lib_list@
smooth_CCLD              = $(F77)
dxmath_SOURCES           = dxmath.c 
dxmath_LDADD             = ../../src/aaa_lib/libapbs.la @lib_list@
dxmath_LDFLAGS           = @LDFLAGS@ @FFLAGS@
dxmath_CCLD              = $(F77)
uhbd_asc2bin_SOURCES     = uhbd_asc2bin.f
uhbd_asc2bin_LDADD       = ../../src/aaa_lib/libapbs.la @lib_list@
uhbd_asc2bin_LDFLAGS     = @LDFLAGS@ @FFLAGS@
uhbd_asc2bin_CCLD        = $(F77)
value_SOURCES            = value.c 
value_LDADD              = ../../src/aaa_lib/libapbs.la @lib_list@
value_LDFLAGS            = @LDFLAGS@ @FFLAGS@
value_CCLD               = $(F77)
benchmark_SOURCES = benchmark.c 
benchmark_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@
similarity_SOURCES = similarity.c 
similarity_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@
analysis_SOURCES = analysis.c 
analysis_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@
dx2mol_SOURCES = dx2mol.c 
dx2mol_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@

