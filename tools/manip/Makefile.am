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

noinst_PROGRAMS  = acc coulomb
acc_SOURCES      = acc.c 
acc_LDADD        = ../../src/aaa_lib/libapbs.la @lib_list@
coulomb_SOURCES  = coulomb.c 
coulomb_LDADD    = ../../src/aaa_lib/libapbs.la @lib_list@
INCLUDES         = @inc_list@
CFLAGS           = @CFLAGS@
