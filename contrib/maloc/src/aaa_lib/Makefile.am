### ##########################################################################
### MALOC = < Minimal Abstraction Layer for Object-oriented C >
### Copyright (C) 1994--2000  Michael Holst
###
### This program is free software; you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation; either version 2 of the License, or (at your
### option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
### See the GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License along
### with this program; if not, write to the Free Software Foundation, Inc.,
### 675 Mass Ave, Cambridge, MA 02139, USA.
###
### rcsid="$Id: Makefile.am,v 1.48 2006/06/03 07:22:28 mholst Exp $"
### ##########################################################################

## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will use the following:
##
##     COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
##     LINK = $(CC) $(CFLAGS) $(LDFLAGS) -o $@
##
## Author:  Michael Holst
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .

BASE_LIBS  = ../base/libbase.la
VSYS_LIBS  = ../vsys/libvsys.la
VSH_LIBS   = ../vsh/libvsh.la
PSH_LIBS   = ../psh/libpsh.la

EFENCE_LIBS = ../efence/libefence.la
if BUILD_EFENCE
LIBEFENCE = libefence.la
else
LIBEFENCE =
endif

libdir = ${prefix}/lib
lib_LTLIBRARIES = libmaloc.la ${LIBEFENCE}

libmaloc_la_SOURCES =
libmaloc_la_LIBADD  = ${BASE_LIBS} ${VSYS_LIBS} ${VSH_LIBS} ${PSH_LIBS}
libmaloc_la_LDFLAGS = -version-info ${FETK_VERSION}

libefence_la_SOURCES =
libefence_la_LIBADD  = ${EFENCE_LIBS}
libefence_la_LDFLAGS = -version-info ${FETK_VERSION}

