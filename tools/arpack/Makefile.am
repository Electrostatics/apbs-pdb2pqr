## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .c .f .o
CCLD = $(F77)

noinst_PROGRAMS  = driver
driver_SOURCES   = driver.c readhb.f ardriv.f
driver_LDADD     = ../../src/aaa_lib/libapbs.la @arpack_lib@ @apbs_build_libraries@
INCLUDES         = @inc_list@
AM_CFLAGS           = @CFLAGS@
AM_LDFLAGS           = @LDFLAGS@ 
