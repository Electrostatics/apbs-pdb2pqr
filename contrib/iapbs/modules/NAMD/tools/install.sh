#!/bin/sh
#
# a script to install iAPBS/NAMD with Intel compilers 11.1
# on Linux x86_64
#
# rok 111809

# create a build directory APBS_PREFIX and download APBS, iAPBS 
# and NAMD source code tarballs there.
export CC=icc
export F77=ifort
export APBS_PREFIX=`pwd`
export APBS_LIB=${APBS_PREFIX}/lib

# unpack the source
tar xzf apbs-1.2.0.tgz
tar xzf iapbs-1.0.0.tgz
tar xzf NAMD_2.7b2_Source.tar.gz

# apbs build
cd $APBS_PREFIX/apbs
export APBS_SRC=`pwd`
./configure --prefix=${APBS_PREFIX} --enable-abps-quiet --disable-openmp
make && make install

cd $APBS_PREFIX/include
ln -s $APBS_SRC/contrib/include/maloc
cd $APBS_PREFIX/lib
cp $APBS_SRC/contrib/lib/libmaloc.a .
cp $APBS_SRC/contrib/blas/.libs/libapbsblas.a .

# iapbs build
cd $APBS_PREFIX/iapbs
./setup.sh
./configure --prefix=${APBS_PREFIX} --disable-openmp
make && make install


# fftw and charm build
cd $APBS_PREFIX/NAMD_2.7b1_Source
export NAMD_SRC=`pwd`

tar xf charm-6.1.3.tar
cd charm-6.1.3

./build charm++ net-linux-x86_64 tcp icc --no-build-shared -O -DCMK_OPTIMIZE=1

cd $NAMD_SRC
wget -O - \
 http://www.ks.uiuc.edu/Research/namd/libraries/fftw-linux-x86_64.tar.gz \
 | tar xzf -
mv linux-x86_64 fftw

wget -O - \
 http://www.ks.uiuc.edu/Research/namd/libraries/tcl-linux-x86_64.tar.gz \
 | tar xzf -
mv linux-x86_64 tcl

# namd2 build
cd $NAMD_SRC
# for NAMD 2.7b2 source code
patch -p1  < $APBS_PREFIX/iapbs/modules/NAMD/patches/namd2-apbs-2.7b2.patch
# if using NAMD from CVS use this patch instead:
#patch -p0 < $APBS_PREFIX/iapbs/modules/NAMD/patches/namd2-apbs-20091117.patch
cd src
ln -s $APBS_PREFIX/iapbs/modules/NAMD/GlobalMasterAPBS.* .
cd ../arch
ln -s $APBS_PREFIX/iapbs/modules/NAMD/*.apbs .

cd $NAMD_SRC
cat >Make.charm <<END
CHARMBASE = $NAMD_SRC/charm-6.1.3
END

cd arch
arch=Linux-x86_64
cat >${arch}.fftw <<END
FFTDIR=${NAMD_SRC}/fftw
FFTINCL=-I\$(FFTDIR)/include
FFTLIB=-L\$(FFTDIR)/lib -lsrfftw -lsfftw
FFTFLAGS=-DNAMD_FFTW
FFT=\$(FFTINCL) \$(FFTFLAGS)
END

cat > ${arch}.tcl <<END
TCLDIR=${NAMD_SRC}/tcl
TCLINCL=-I\$(TCLDIR)/include
TCLLIB=-L\$(TCLDIR)/lib -ltcl8.3 -ldl
TCLFLAGS=-DNAMD_TCL
TCL=\$(TCLINCL) \$(TCLFLAGS)
END

cd $NAMD_SRC

./config apbs ${arch}-icc --charm-arch net-linux-x86_64-tcp-icc
cd ${arch}-icc
make -j4


# to test
./namd2 src/alanin
./charmrun ++local +p2 ./namd2 src/alanin

cd $APBS_PREFIX/iapbs/modules/NAMD/test/dipeptide
$NAMD_SRC/${arch}-icc/namd2 apbs-solvation.conf| tee apbs-solvation.out
$NAMD_SRC/${arch}-icc/namd2 apbs-md.conf| tee apbs-md.out

