## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .c .f .o

always_built_SUBDIRS = @apbs_do_blas@ @apbs_do_maloc@ pmgZ aqua opal @zlib@

SUBDIRS = $(always_built_SUBDIRS) .
