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
if [ $NCF -eq 1 ]; then
  cd $NCDIR
  mkdir -p netcdf-f
  cd netcdf-f
#  wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.2/netcdf-fortran-4.6.2.tar.gz
  tar xzf netcdf-fortran-4.6.2.tar.gz
  cd netcdf-fortran-4.6.2
  CPPFLAGS="-I$NCCDIR/netcdf-c-4.9.3/include"
  LDFLAGS="-L$NCCDIR/netcdf-c-4.9.3/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lnetcdf"
  FC=pgf90
  CC=gcc
  ./configure --prefix=$NCFDIR

  make -j4 check
  make install
fi

