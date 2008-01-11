## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .cc .c .f .o
CCLD = $(F77)

build_libraries = @apbs_build_libraries@

bindir = ${prefix}/bin
bin_PROGRAMS     = apbs
apbs_SOURCES     = main.c 
apbs_LDADD       = ${build_libraries}

libdir = ${prefix}/lib
lib_LTLIBRARIES = libapbsmainroutines.la
libapbsmainroutines_la_SOURCES = routines.c
libapbsmainroutines_la_LDFLAGS = -version-info ${APBS_VERSION}

INCLUDES = @inc_list@
AM_CFLAGS= @wall_pedantic@
AM_LDFLAGS       = @LDFLAGS@ @FFLAGS@
 
includedir = ${prefix}/include/apbs
include_HEADERS = routines.h
