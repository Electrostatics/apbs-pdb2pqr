#! /bin/sh

##############################################################################
# File:    bootstrap
#
# Purpose: Auto-generate the GNU configure script.
#          (WITH LIBTOOL AND AUTOHEADER.)
#
# Author:  Michael Holst
##############################################################################

rm -rf config.cache autom4te.cache

aclocal --verbose 
automake --verbose --gnu --add-missing --copy
autoconf --verbose 
autoheader --verbose 

if [ -x libtoolize ]; then
  libtoolize --automake --copy --force
fi

rm -rf config.cache autom4te.cache

