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
THISLIB = fem

noinst_LTLIBRARIES = libapbsfem.la
libapbsfem_la_SOURCES  = vfetk.c vcsm.c vpee.c dummy.c
INCLUDES           = -I${top_srcdir}/src/generic -I${top_srcdir}/src/mg @inc_list@
AM_CFLAGS          = @wall_pedantic_ansi@ @CFLAGS@


all-local:
	-cp -p ${top_srcdir}/src/${THISLIB}/apbs/*.h ${top_srcdir}/src/aaa_inc/apbs/.

