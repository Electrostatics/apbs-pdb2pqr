
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
## Author:  Michael Holst
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .cc .c .f .o
INCLUDES           = @inc_list@
LDFLAGS            = @LDFLAGS@ @FLIBS@

noinst_PROGRAMS    = mgmesh 
mgmesh_SOURCES     = mgmesh.c mypde.f
mgmesh_LDADD       = ../../src/aaa_lib/libapbs.la @lib_list@

