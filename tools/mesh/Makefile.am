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
uhbd_asc2bin_LDADD       = ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
uhbd_asc2bin_LDFLAGS     = @LDFLAGS@ @FFLAGS@
uhbd_asc2bin_CCLD        = $(F77)
else
UHBD =
endif
noinst_PROGRAMS          = mgmesh dxmath mergedx mergedx2 $(UHBD) smooth value dx2mol dx2uhbd benchmark similarity analysis multivalue del2dx tensor2dx

mgmesh_SOURCES				= mgmesh.c 
mgmesh_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
mgmesh_LDFLAGS				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
mgmesh_CCLD					= $(CC)
smooth_SOURCES				= smooth.c 
smooth_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
smooth_LDFLAGS				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
smooth_CCLD					= $(F77)
dxmath_SOURCES				= dxmath.c 
dxmath_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
dxmath_LDFLAGS				= @LDFLAGS@ 
dxmath_CCLD					= $(F77)
mergedx_SOURCES				= mergedx.c
mergedx_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
mergedx_LDFLAGS				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
mergedx_CCLD				= $(CC)
mergedx2_SOURCES			= mergedx2.c
mergedx2_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
mergedx2_LDFLAGS			= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
mergedx2_CCLD				= $(CC)
value_SOURCES				= value.c 
value_LDADD					= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
value_LDFLAGS				= @LDFLAGS@ 
value_CCLD					= $(F77)
multivalue_SOURCES			= multivalue.c 
multivalue_LDADD			= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
multivalue_LDFLAGS			= @LDFLAGS@ 
multivalue_CCLD				= $(F77)
benchmark_SOURCES			= benchmark.c 
benchmark_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
similarity_SOURCES			= similarity.c 
similarity_LDADD			= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
analysis_SOURCES			= analysis.c 
analysis_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
dx2mol_SOURCES				= dx2mol.c 
dx2mol_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
dx2uhbd_SOURCES				= dx2uhbd.c 
dx2uhbd_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
del2dx_SOURCES				= del2dx.c 
del2dx_LDADD				= ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
tensor2dx_SOURCES                       = tensor2dx.c
tensor2dx_LDADD                         = ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
