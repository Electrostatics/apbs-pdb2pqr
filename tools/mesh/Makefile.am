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

noinst_PROGRAMS          = mgmesh dxmath uhbd_asc2bin smooth
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
dxmath_LDFLAGS           = @LDFLAGS@ @FLIBS@
dxmath_CCLD              = $(F77)
uhbd_asc2bin_SOURCES     = uhbd_asc2bin.f
uhbd_asc2bin_LDADD       = ../../src/aaa_lib/libapbs.la @lib_list@
uhbd_asc2bin_LDFLAGS     = @LDFLAGS@
uhbd_asc2bin_CCLD        = $(F77)
