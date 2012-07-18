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
THISLIB = pmgc

noinst_LTLIBRARIES = libapbspmgc.la
libapbspmgc_la_SOURCES  = buildAd.c buildBd.c buildGd.c buildPd.c cgd.c gsd.c \
			  matvecd.c mgcsd.c mgdrvd.c mgsubd.c mikpckd.c \
			  mlinpckd.c mypdec.c newtond.c newdrvd.c powerd.c \
			  smoothd.c mgfasd.c

INCLUDES           = -I${top_srcdir}/src/generic @inc_list@
AM_CFLAGS          = @wall_pedantic_ansi@ @CFLAGS@


all-local:
	-cp -p ${top_srcdir}/src/${THISLIB}/apbs/*.h ${top_srcdir}/src/aaa_inc/apbs/.

