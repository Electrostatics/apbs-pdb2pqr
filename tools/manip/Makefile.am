## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .c .f .o
CCLD = $(CC)

noinst_PROGRAMS  = coulomb born
coulomb_SOURCES  = coulomb.c 
coulomb_LDADD    = ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
born_SOURCES     = born.c 
born_LDADD       = ../../src/aaa_lib/libapbs.la @apbs_build_libraries@
INCLUDES         = @inc_list@
AM_CFLAGS        = @CFLAGS@
AM_LDFLAGS       = @LDFLAGS@
