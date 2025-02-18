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

```
Program Name:assl2
Data A is <time series A text file>
Data B is <time series B text file>

Number of the data set A = <N>
Number of the data set B = <N>
Options
Mode                        =  <0 or 1>
Initial step of data A      =  <n1>
Initial step of data B      =  <n2>
Final   step of data A      =  <N1>
Final   step of data B      =  <N2>
egrees of freedom           =  <M>
alue of sparsity parameter  =  <r>          

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
...
The Sparse structure of data set B
  1   2
  1   3
  1   4
  1   5
  3   4
  3   5
  3  13
  4   5
  4   6
  4   7
...
The anomaly of each dimension
  #      a->b     b->a     max)
  1    0.000    0.000    0.000
  2    0.000    0.000    0.000
  3    0.000    0.000    0.000
  4    0.000    0.000    0.000
  5    0.001    0.001    0.001
  6    0.001    0.001    0.001
  7    0.000    0.000    0.000
  8    0.000    0.000    0.000
  9    0.000    0.000    0.000
 10    0.000    0.000    0.000
```

This output consists of several blocks, each presenting different information:

The first block displays the input file.
The second block shows the contents of the options used.
The third block presents the convergence of the computation.
The block "The sparse structure of data set A" illustrates the sparse correlation of data set A.
The paired numbers represent the indices of degrees of freedom that exhibit correlation. Similarly, the block "The sparse structure of data set B" conveys the same meaning for data set B.

Finally, the block "The anomaly of each dimension" identifies the degrees of freedom that exhibit anomalies by comparing data set A and data set B.
The first column represents the index of the degree of freedom, and when the value in the fourth column is nonzero, that degree of freedom is considered to have an anomaly.

