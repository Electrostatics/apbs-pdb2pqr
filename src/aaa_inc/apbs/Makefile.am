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

BOTH_HEADS  = apolparm.h femparm.h mgparm.h nosh.h pbeparm.h vacc.h valist.h \
              vatom.h vcap.h vgreen.h vhal.h vpbe.h vunit.h vstring.h \
              vparam.h vclist.h
PMG_HEADS   = vpmg.h vpmgp.h vgrid.h vmgrid.h vopot.h
FETK_HEADS  = vfetk.h vcsm.h vpee.h

APBS_HEADS = apbs.h

OTHER_HEADS = ${BOTH_HEADS} ${PMG_HEADS} ${FETK_HEADS}

includedir = ${prefix}/include/apbs
include_HEADERS = ${APBS_HEADS} ${OTHER_HEADS}

distclean-local: clean-local

clean-local:
	-rm -f ${OTHER_HEADS}

