### ##########################################################################
### MC = < Manifold Code >
### Copyright (C) 1993  Michael Holst
### Copyright (C) 1999  Contributing authors listed in the code documentation
###                     retain the copyrights to their code contributions.
###
### The program MC is copyrighted software, and there are restrictions on its
### use, modification, and distribution; please refer to the document COPYING
### which accompanies this software.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
### See the document COPYING for more details.
###
### You should have received a copy of the COPYING file with this program;
### if not, write to the author Michael Holst at:  mholst@math.ucsd.edu
###
### rcsid="$Id$"
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
## Author:  Nathan Baker and Michael Holst
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .c .f .o

always_built_SUBDIRS = vatom vhash vcsm vpbe

SUBDIRS = $(always_built_SUBDIRS) @sav@
INCLUDES = @mcdir@/inc
