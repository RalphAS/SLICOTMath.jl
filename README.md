# SLICOTMath
[![GitHub CI Build Status](https://github.com/RalphAS/SLICOTMath.jl/workflows/CI/badge.svg)](https://github.com/RalphAS/SLICOTMath.jl/actions)
[![Coverage Status](http://codecov.io/github/RalphAS/SLICOTMath.jl/coverage.svg?branch=master)](http://codecov.io/github/RalphAS/SLICOTMath.jl?branch=master)

# Introduction
This package provides Julia wrappers for some routines in the SLICOT library,
which contains Fortran implementations of control theory algorithms.

SLICOT is hosted at https://github.com//SLICOT/SLICOT-Reference.git

The group of mathematical routines (those starting with "M") are covered here,
except for a few which are awkward to call from Julia (there are preferable implementations
of those algorithms elsewhere in the Julia ecosystem).

For the most part the Julia API resembles that of LAPACK routines in the LinearAlgebra
standard library. Thus array dimension arguments are (almost always) suppressed and
workspace is usually allocated as needed within the wrappers. One difference is that
result arrays are only passed as arguments to be mutated, and are not among the explicitly
returned objects.

# Status: Caveat Emptor
This is work in progress. So far all of the wrappers and tests were mechanically generated
and edited to repair oversights.

For substantive documentation the user is referred to the HTML documents
in [the SLICOT repository](https://github.com//SLICOT/SLICOT-Reference.git)

Users will also need to refer to the method signatures. Argument names are lower case
versions of those in the wrapped routines.

## Tests

The test routines generally verify that the API is consistent. They were mechanically
translated from the SLICOT-Reference examples, so they include extra computations (beyond
merely invoking the wrappers).  The output of most has been visually compared to the
analogous Fortran results, almost always agreeing (or differing in ways consistent with
floating-point and BLAS/LAPACK version variations).  That comparison is *not* currently
part of the package test framework. The CI tests as such rarely check the accuracy of
algorithms, although residuals are occasionally reported.

Interfaces to routines which are not called by the reference example codes have been
quickly reviewed, but very few have been tested by the author.


# Related work
Many of the higher-level control-theoretic methods in SLICOT have been translated to
Julia in [packages written by Andreas Varga](https://github.com//AndreasVarga).

No endorsement or promotion of SLICOTMath.jl by the authors of SLICOT
is implied.
