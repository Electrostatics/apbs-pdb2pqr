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
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .
THISLIB = generic

noinst_LTLIBRARIES = libapbsgen.la
libapbsgen_la_SOURCES  = nosh.c mgparm.c femparm.c vacc.c valist.c vatom.c \
                         vgreen.c vpbe.c vcap.c
INCLUDES           = -I../fem -I../mg @inc_list@ 
AM_CFLAGS          = @wall_pedantic_ansi@
LDFLAGS          = @LDFLAGS@ @FLIBS@


all-local:
	-cp -p ${top_srcdir}/src/${THISLIB}/apbs/*.h ../aaa_inc/apbs/.

