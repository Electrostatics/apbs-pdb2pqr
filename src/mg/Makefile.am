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
libapbsmg_la_SOURCES  = vpmg.c vpmgp.c buildAd.f buildBd.f buildGd.f \
  buildPd.f cgd.f cgmgd.f cgmgdrvd.f dlphsubd.f gsd.f mainsubd.f matvecd.f \
  mblasd.f mgcsd.f mgdrvd.f mgfasd.f mgsubd.f mikpckd.f mlinpckd.f ncgd.f \
  ncgdrvd.f newdrvd.f newtond.f ngsd.f ngsdrvd.f ninterpd.f nrchd.f \
  nrchdrvd.f nsmoothd.f nsord.f nsordrvd.f nwjd.f nwjdrvd.f powerd.f rchd.f \
  smoothd.f sord.f wjd.f mypde.f vpmgf.f
INCLUDES           = -I../generic -I../fem @inc_list@
AM_CFLAGS          = @wall_pedantic_ansi@
CFLAGS             = @CFLAGS@
LDFLAGS          = @LDFLAGS@ @FLIBS@


all-local:
	-cp -p ${top_srcdir}/src/${THISLIB}/apbs/*.h ../aaa_inc/apbs/.

