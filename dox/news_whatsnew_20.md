What's new in release 2.0 {#news_whatsnew_20}
=========================

\page news_whatsnew_20 What's new in release 2.0?

**Fast multipole methods**

The main new feature of the 2.0 release is the inclusion of the generic fast multipole extension.
The key characteristics of the extension are:

- Generic implementation of the fast multipole method (FMM) that can be specialized for solving problems of various types
- Cluster tree building based on specific cluster division criteria
- FMM operator manipulations including arithmetic expressions
- Transparent interface of the generalized FMM matrix class applicable in iterative linear and eigensolvers
- Parallelization using the OpenMP standard, BFS and DFS traversing strategies, built-in timing analysis module
- Contains specialized algorithms for different fundamental solutions, as well as general black-box method for smooth kernels
- Accessible both in standalone and Matlab-based applications

**Nearly singular integrals**

New methods are integrated to handle nearly singular integrals arising in various boundary element formulations.

- Can be specialized for various kernel and field types
- Supports compile-time algorithm selection
- Contains the implementation of the adaptive quadrature method proposed by Telles.

**Applications**

The new release contains a collection of prepared applicaions that can applied to solve problems in different fields of engineering.
The prepared applications also serve as guidelines in developing new, user-customized solvers.
Some applications are listed below.

- Elastodynamics and elastostatics: 3D conventional BEM
- Helmholtz solvers: 2D conventional BEM and wideband FMBEM, 3D conventional BEM and high frequency FMBEM, Burton & Miller formulation for mitigating the problem of fictitious eigenfrequencies
- Laplace solvers: 2D conventional BEM and FMBEM, 3D conventional BEM and FMBEM
- Stochastic eigendecomposition: Generic Black-box FMM (arbitrary space and field dimensions)

For a detailed list of collected appliactions please visit the \ref applications page.

**Benchmark cases for linear acoustics**

Release 2.0 contains some selected benchmark problems in the field of linear acoustics.


