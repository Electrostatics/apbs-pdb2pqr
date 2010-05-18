## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  David Gohara
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .c .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .
THISLIB = z

noinst_LTLIBRARIES = libz.la
libz_la_SOURCES = adler32.c compress.c crc32.c deflate.c gzclose.c \
	gzlib.c gzread.c gzwrite.c infback.c inffast.c inflate.c inftrees.c trees.c zutil.c
	
AM_FLAGS			  = @CFLAGS@
INCLUDES           = -I${top_srcdir}/contrib/zlib