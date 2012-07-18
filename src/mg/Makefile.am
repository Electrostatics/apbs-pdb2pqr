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
THISLIB = mg

noinst_LTLIBRARIES = libapbsmg.la
libapbsmg_la_SOURCES = vgrid.c vopot.c vpmg.c vpmgp.c

INCLUDES = -I${top_srcdir}/src/generic -I${top_srcdir}/src/pmgc -I${top_srcdir}/src/fem @inc_list@
AM_CFLAGS = @wall_pedantic_ansi@ @CFLAGS@


all-local:
	-cp -p ${top_srcdir}/src/${THISLIB}/apbs/*.h ${top_srcdir}/src/aaa_inc/apbs/.

