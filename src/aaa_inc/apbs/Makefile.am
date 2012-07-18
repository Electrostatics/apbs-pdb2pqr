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

GENERIC_HEADS  = apolparm.h femparm.h mgparm.h nosh.h pbeparm.h vacc.h \
		 valist.h vatom.h vcap.h vgreen.h vhal.h vpbe.h vunit.h \
		 vstring.h vparam.h vclist.h vmatrix.h

PMGC_HEADS     = buildAd.h buildBd.h buildGd.h buildPd.h cgd.h gsd.h matvecd.h \
		 mgcsd.h mgdrvd.h mgfasd.h mgsubd.h mikpckd.h mlinpckd.h \
		 mypdec.h newdrvd.h newtond.h powerd.h smoothd.h

MG_HEADS       = vpmg.h vpmgp.h vgrid.h vmgrid.h vopot.h

FETK_HEADS = vfetk.h vcsm.h vpee.h

APBS_HEADS = apbs.h

OTHER_HEADS = ${GENERIC_HEADS} ${PMGC_HEADS} ${MG_HEADS} ${FETK_HEADS}

includedir = ${prefix}/include/apbs
include_HEADERS = ${APBS_HEADS} ${OTHER_HEADS}

distclean-local: clean-local

clean-local:
	-rm -f ${OTHER_HEADS}

