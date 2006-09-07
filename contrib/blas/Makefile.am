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
THISLIB = apbsblas

lib_LTLIBRARIES = ${prefix}/libapbsblas.la
libapbsblas_la_SOURCES  = mblasd.f
INCLUDES           = @inc_list@
AM_CFLAGS          = @wall_pedantic_ansi@ @CFLAGS@

