### ##########################################################################
# APBS -- Adaptive Poisson-Boltzmann Solver
#
#  Nathan A. Baker (nbaker@wasabi.ucsd.edu)
#  Dept. of Chemistry and Biochemistry
#  Dept. of Mathematics, Scientific Computing Group
#  University of California, San Diego 
#
#  Additional contributing authors listed in the code documentation.
#
# Copyright © 1999. The Regents of the University of California (Regents).
# All Rights Reserved. 
# 
# Permission to use, copy, modify, and distribute this software and its
# documentation for educational, research, and not-for-profit purposes,
# without fee and without a signed licensing agreement, is hereby granted,
# provided that the above copyright notice, this paragraph and the
# following two paragraphs appear in all copies, modifications, and
# distributions.
# 
# IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
# SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
# ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
# REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
# 
# REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
# ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
# TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS. 
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
SUFFIXES = .H .h .C .cc .c .f .o

bin_PROGRAMS     = apbs
apbs_SOURCES     = main.c mypde.f
apbs_LDADD       = -lapbs @lib_list@
INCLUDES         = @inc_list@
LDFLAGS          = @LDFLAGS@ @FLIBS@
