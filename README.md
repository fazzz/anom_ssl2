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

[-mode ]   `if 0 normal standarization will be done, if 1 not`

[-ninia]    `initial step of data A`

[-ninib]    `initial step of data B`

[-nfina]    `final   step of data A`

[-nfinb]    `final   step of data B`

[-Sparsity] `sparsity (from 0 to 1)`

[-ndim]     `dimensionality of data`

## Output

           Program Name:assl2
        Data A is ../dist_r-r_3/r-r_dist_union-open-close.txt
        Data B is ../../../close/analyses_2/dist_r-r_3/r-r_dist_union-open-close.txt

       Number of the data set A = ****
       Number of the data set B = ****
        Options
        Mode                        =            0
        Initial step of data A      =            1
        Initial step of data B      =            1
        Final   step of data A      =        10001
        Final   step of data B      =        10001
        Degrees of freedom          =          912
        Value of sparsity parameter =    0.900

The block coordinate descent method will be started.
  2-th max diff:  0.020183
  3-th max diff:  0.000622
  4-th max diff:  0.000005
  5-th max diff:  0.000000
  6-th max diff:  0.000000
  7-th max diff:  0.000000
The block coordinate descent method is converged til   8-th iteration.
The block coordinate descent method will be started.
  2-th max diff:  0.017372
  3-th max diff:  0.000449
  4-th max diff:  0.000006
  5-th max diff:  0.000000
  6-th max diff:  0.000000
The block coordinate descent method is converged til   7-th iteration.
The Sparse structure of data set A
  1   2
  1   3
  1   4
  1   5
  1   6
  1   9
  2   3
  2   4
  2   5
  2   6
  2   9




