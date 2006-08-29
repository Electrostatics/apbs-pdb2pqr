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
### rcsid="$Id: Makefile.am,v 1.46 2006/06/03 07:22:28 mholst Exp $"
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

VSYS_HEADS   = vsys.h   vio.h vmem.h vnm.h vset.h
VSH_HEADS    = vsh.h
PSH_HEADS    = psh.h   vmp.h vcom.h vmpi.h

if BUILD_EFENCE
EFENCE_HEADS = 
else
EFENCE_HEADS = 
endif

MALOC_HEADS = maloc.h maloc_base.h
OTHER_HEADS = ${VSYS_HEADS} ${VSH_HEADS} ${PSH_HEADS} ${EFENCE_HEADS}

includedir = ${prefix}/include/maloc
include_HEADERS = ${MALOC_HEADS} ${OTHER_HEADS}

distclean-local: clean-local

clean-local:
	-rm -f ${OTHER_HEADS}

