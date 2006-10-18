## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .
THISLIB = generic

noinst_LTLIBRARIES = libapbsgen.la
libapbsgen_la_SOURCES  = nosh.c mgparm.c femparm.c pbeparm.c apolparm.c \
                         vacc.c valist.c vatom.c vgreen.c vpbe.c vcap.c \
                         vclist.c vstring.c vparam.c 
INCLUDES           = -I${top_srcdir}/src/fem -I${top_srcdir}/src/mg -I${top_srcdir}/src/aaa_inc @inc_list@ 
AM_CFLAGS          = @wall_pedantic_ansi@ @CFLAGS@

all-local:
	-cp -p ${top_srcdir}/src/${THISLIB}/apbs/*.h ${top_srcdir}/src/aaa_inc/apbs/.
