## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS = apbs
SUBDIRS = $(always_built_SUBDIRS) .


includedir = ${prefix}/include/apbs
include_HEADERS = apbscfg.h
