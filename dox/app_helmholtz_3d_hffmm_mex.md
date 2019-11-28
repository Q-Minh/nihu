3D Helmholtz FMBEM with Mex interface {#app_helmholtz_3d_hffmm_mex}
=====================================

\page app_helmholtz_3d_hffmm_mex 3D Helmholtz FMBEM with Mex interface

[TOC]

Description {#app_helmholtz_3d_hffmm_mex_desc}
===========

This application implements a fast multipole solver for the Helmholtz equation in 3D with the following features:

- Collocational formulation
- Homogeneous meshes with quadrilateral elements
- High frequency multilevel fast multipole method
- Burton-Miller technique for mitigating fictitious eigenfrequencies
- Applicable for Neumann and Dirichlet problems
- MEX interface for interacting with Matlab - Compatible with Matlab's built-in iterative 
solvers
 
The application is used for assembling the single and double layer potential matrices (SLP and DLP matrices) and then evaluate fast matrix-vector products with these matrices.
The same matrix-vector products can be called from Matlab's built-in iterative solvers, such as BiCGSTAB or GMRES.

Usage {#app_helmholtz_3d_hffmm_mex_usage}
=====

This MEX application can be executed in various steps.
Each step is represented by a command string that is passed to the MEX function as its first parameter.
The corresponding MEX function has an inner state, therefore initialization and cleanup are necessary.
The typical sequence of the steps and the corresponding commands is as follows.

1. `'init'` -- Initializes the data structure and the storage object.
2. `'set'` -- Sets up the parameters for the computation, such as the wave number, or the multipole expansion accuracy.
3. `'mesh'` -- Defines the surface and field meshes. The field mesh is not necessarily given, in this case the surface problem is solved with singular integral evaluation and the Burton-Miller technique. If the field mesh is also given, the standard collocational method is applied for computing the radiated field.
4. `'tree'` -- Build the cluster tree and initialize the interaction lists of the clusters.
5. `'matrix'` -- Assemble the FMM matrices for subsequent fast matrix-vector products. This step includes the precomputation of the discretized operators.
6. `'mvp_dlp'` and `'mvp_slp'` Evaluate fast matrix-vector product by the SLP or the DLP matrix. These steps are called in each iteration of an iterative solver. The same steps are called once for computing the radiated field.
7. `'cleanup'` -- Free the memory allocated for the internal storage and reset to the initial state.

Arriving at a solution of a complete Helmholtz problem may require solving the surface problem at first, and then computing the radiated field.
In this case, the whole sequence of steps is repeated twice, first the appilaction is used in surface mode, and then in field mode.
The example in the following section also demonstrates this case.

Example {#app_helmholtz_3d_hffmm_mex_app}
=======

This example demonstrates the usage of the application through a transparent problem.
The field of a point source is computed by prescribing the normal derivative field of the source on a surface mesh as a Neumann boundary condition.
The solution is performed in two steps: first the field on the surface is attained by solving the BEM system of equations using the fast multipole method, and then, the field in the field points is evaluated by means of fast matrix-vector products.
The numerical solution of the transparent problem is finally compared to the analytical solution.


