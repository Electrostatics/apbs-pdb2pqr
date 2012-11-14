#! /bin/sh
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
### rcsid="$Id: setup.sh 1125 2007-09-06 20:52:30Z sdg0919 $"
### ##########################################################################

##############################################################################
# File:    bootstrap
#
# Purpose: Auto-generate the GNU configure script.
#          (WITH LIBTOOL AND AUTOHEADER.)
#
# Author:  Michael Holst
##############################################################################

# rm -rf aclocal.m4  autom4te.cache/ config config.status config.log configure libtool m4/ Makefile Makefile.in
# autoreconf --warnings=all --force --verbose --install

if [ ! -x configure ] ; then
 echo "Creating configure script ..."
 rm -rf config.cache autom4te.cache
 mkdir config
 aclocal \
 && automake --gnu --add-missing --copy \
 && autoconf \
 && libtoolize --automake --copy --force
 rm -rf config.cache autom4te.cache
fi

#aclocal --force
#mkdir config
#/usr/bin/autoconf --force --warnings=all
#automake --add-missing --copy --force-missing --warnings=all

# setup libs
echo "Setting up library dependencies ..."
p=`pwd`
export APBS_SRC=`echo $p| sed 's%/contrib/iapbs%%'`
export APBS_PREFIX=`grep "^prefix =" ../../Makefile | awk '{print $3}'`

cp -a $APBS_SRC/contrib/include/maloc $APBS_PREFIX/include
cp $APBS_SRC/contrib/lib/libmaloc.a $APBS_PREFIX/lib
cp $APBS_SRC/contrib/blas/.libs/libapbsblas.a $APBS_PREFIX/lib
cp $APBS_SRC/contrib/zlib/.libs/libz.a $APBS_PREFIX/lib

#echo "APBS_SRC=$APBS_SRC"
#echo "APBS_PREFIX=$APBS_PREFIX"

#echo "Running configure and building iAPBS ..."
#./configure --prefix=${APBS_PREFIX} --disable-openmp \
#&& make && make install

