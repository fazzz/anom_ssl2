# Anomaly detection

This is a tool based on Sparse Structure Learning (Ide et. al. 2008) for Givin Two Data sets of Time Serie

<!--
## The algorithms
-->

## Requirement

LAPACK (https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz)

<!--
``
wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz
tar -zxvf v3.10.0.tar.gz
cd lapack-3.10.0
cp make.inc.example make.inc
emacs make.inc
make -j 4
cp lapack-3.10.0/liblapack.a $LAPACK/
cp librefblas.a $LAPACK/libblas.a
cp lapack-3.10.0/libtmglib.a $LAPACK/
``
-->

## Installation

``
./configure LDFLAGS=-$LAPACK/lib
make && make install
``

## How to use

``
assl2 -Sparsity *r* -ndim *ndist* *time_series_A* *time_series_B* > *out*
``

## Options

The followings are the options for `assl2`.

  -mode     :  if 0 normal standarization will be done, if 1 not
  -ninia    :  initial step of data A
  -ninib    :  initial step of data B
  -nfina    :  final   step of data A
  -nfinb    :  final   step of data B
  -Sparsity :  sparsity
  -ndim     :  dimensionality of data

