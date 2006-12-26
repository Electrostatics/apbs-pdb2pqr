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

if UHBD_COMPILE
UHBD                     = uhbd_asc2bin
uhbd_asc2bin_SOURCES     = uhbd_asc2bin.f
uhbd_asc2bin_LDADD       = ../../src/aaa_lib/libapbs.la @lib_list@
uhbd_asc2bin_LDFLAGS     = @LDFLAGS@ @FFLAGS@
uhbd_asc2bin_CCLD        = $(F77)
else
UHBD =
endif
noinst_PROGRAMS          = mgmesh dxmath mergedx $(UHBD) smooth value dx2mol dx2uhbd benchmark similarity analysis multivalue
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
mergedx_SOURCES          = mergedx.c
mergedx_LDADD            = ../../src/aaa_lib/libapbs.la @lib_list@
mergedx_LDFLAGS          = ../../src/aaa_lib/libapbs.la @lib_list@
mergedx_CCLD             = $(CC)
value_SOURCES            = value.c 
value_LDADD              = ../../src/aaa_lib/libapbs.la @lib_list@
value_LDFLAGS            = @LDFLAGS@ @FFLAGS@
value_CCLD               = $(F77)
multivalue_SOURCES       = multivalue.c 
multivalue_LDADD         = ../../src/aaa_lib/libapbs.la @lib_list@
multivalue_LDFLAGS       = @LDFLAGS@ @FFLAGS@
multivalue_CCLD          = $(F77)
benchmark_SOURCES = benchmark.c 
benchmark_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@
similarity_SOURCES = similarity.c 
similarity_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@
analysis_SOURCES = analysis.c 
analysis_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@
dx2mol_SOURCES = dx2mol.c 
dx2mol_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@
dx2uhbd_SOURCES = dx2uhbd.c 
dx2uhbd_LDADD = ../../src/aaa_lib/libapbs.la @lib_list@
