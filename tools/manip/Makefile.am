
## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will use the following:
##
##     COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
##     LINK = $(CC) $(CFLAGS) $(LDFLAGS) -o $@
##
## Author:  Nathan Baker and Michael Holst
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .c .f .o

noinst_PROGRAMS  = acc coulomb
acc_SOURCES      = acc.c 
acc_LDADD        = ../../src/aaa_lib/libapbs.la @lib_list@
coulomb_SOURCES  = coulomb.c 
coulomb_LDADD    = ../../src/aaa_lib/libapbs.la @lib_list@
INCLUDES         = @inc_list@
CFLAGS           = @CFLAGS@
