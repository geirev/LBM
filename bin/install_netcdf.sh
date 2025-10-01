#!/bin/bash

HOME=/home/AD.NORCERESEARCH.NO/geev/
NCDIR=$HOME/.netcdf

ZLIBDIR=$NCDIR/zlib
HDF5DIR=$NCDIR/hdf5
NCCDIR=$NCDIR/netcdf-c
NCFDIR=$NCDIR/netcdf-f

ZLIB=1
HDF5=1
NCC=1
NCF=1

# zlib
if [ $ZLIB -eq 1 ]; then
  cd $NCDIR
  mkdir -p zlib
  cd zlib
  wget https://zlib.net/zlib-1.3.1.tar.gz
  tar xzf zlib-1.3.1.tar.gz
  cd zlib-1.3.1
  ./configure --prefix=$ZLIBDIR
  make -j4 check
  make install
fi

# hdf5
if [ $HDF5 -eq 1 ]; then
  cd $NCDIR
  mkdir -p hdf5
  cd hdf5
  wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz
  tar xzf hdf5-1.14.6.tar.gz
  cd hdf5-1.14.6
  ./configure --with-zlib=$ZLIBDIR --prefix=$HDF5DIR --enable-hl
  make -j4 check
  make install
fi

# netcdf-c
if [ $NCC -eq 1 ]; then
  cd $NCDIR
  mkdir -p netcdf-c
  cd netcdf-c
  wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.3/netcdf-c-4.9.3.tar.gz
  tar xzf netcdf-c-4.9.3.tar.gz
  cd netcdf-c-4.9.3
  CPPFLAGS="-I$HDF5DIR/include -I$ZLIBDIR/include" LDFLAGS="-L$HDF5DIR/lib -L$ZLIBDIR/lib" ./configure --prefix=$NCCDIR
  make -j4 check
  make install
fi

# netcdf-f
if [ $NCF -eq 1 ]; then
  cd $NCDIR
  mkdir -p netcdf-f
  cd netcdf-f
  wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.2/netcdf-fortran-4.6.2.tar.gz
  tar xzf netcdf-fortran-4.6.2.tar.gz
  cd netcdf-fortran-4.6.2
  CPPFLAGS="-I$NCCDIR/include" 
  LDFLAGS="-L$NCCDIR/lib"
  FC=pgf90 
  CC=gcc
  ./configure --prefix=$NCFDIR
  make -j4 check
  make install
fi


## netcdf-f
#if [ $NCF -eq 1 ]; then
#  cd $NCDIR
#  mkdir -p netcdf-f
#  cd netcdf-f
#  wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.2/netcdf-fortran-4.6.2.tar.gz
#  tar xzf netcdf-fortran-4.6.2.tar.gz
#  cd netcdf-fortran-4.6.2
#  CPPFLAGS="-I$NCCDIR/netcdf-c-4.9.3/include"
#  LDFLAGS="-L$NCCDIR/netcdf-c-4.9.3/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lnetcdf"
#  FC=pgf90
#  CC=gcc
#  ./configure --prefix=$NCFDIR
#
#  make -j4 check
#  make install
#fi
