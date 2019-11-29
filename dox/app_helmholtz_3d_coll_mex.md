3D Helmholtz conventional collocational BEM with MEX interface {#app_helmholtz_3d_coll_mex}
==============================================================

\page app_helmholtz_3d_coll_mex 3D Helmholtz conventional collocational BEM with MEX interface

[TOC]

Description {#app_helmholtz_3d_coll_mex_desc}
===========

This application implements the BEM matrix assembly for 3D Helmholtz problems.
The main features of the application are the following:

- Heterogeneous meshes with triangular or quadrilateral elements
- Collocational formalism
- MEX interface for calling the application from Matlab
- Matrices assembled directly into memory allocated by Matlab
 
The assembled surface and field BEM matrices are utilized for solving the BEM system of equations and computing the scattered field at the field points.

Note: this application does not implement any techniques for mitigating fictitious eigenfrequencies.

Usage {#app_helmholtz_3d_coll_mex_usage}
=====

The application can be called from Matlab as

	[Ls, Ms, Lf, Mf] = helmholtz_3d_coll_bem(s_nodes, s_elems, f_nodes, f_elems, k);

Input parameters:

- `s_nodes` -- Real matrix of type \f$ N \times 3 \f$ containing surface node locations (\f$ N \f$ denoting the number of surface nodes)
- `s_elems` -- Real matrix of type \f$ E \times 9 \f$ containing surface element connectivity (\f$ E \f$ denoting the number of elements on the surface). The first column contains the type identifier of each element.
- `f_nodes` -- Real matrix of type \f$ M \times 3 \f$ containing Field mesh node locations, similar to `s_nodes`
- `f_elems` -- Real matrix of type \f$ F \times 9 \f$ containing field mesh element connectivity, similar to `s_elems`
- `k` -- Real scalar wave number

Output:

- `Ls` -- Single layer potential surface matrix (\f$ E \times E \f$ full complex valued matrix)
- `Ms` -- Double layer potential surface matrix (\f$ E \times E \f$ full complex valued matrix)
- `Lf` -- Single layer potential field matrix (\f$ F \times F \f$ full complex valued matrix)
- `Mf` -- Double layer potential field matrix (\f$ F \times F \f$ full complex valued matrix)

For solving the Neumann problem on the surface one can use 

	ps = Ms \ (Ls * qs);

Similarly, the Dirichlet problem is solved as

	qs = Ls \ (Ms * ps);

with `ps` and `qs` denoting the \f$ E \times 1 \f$ complex valued column vectors containing the field and its normal derivative in the surface collocation points, respectively.

Finally, using the vectors `ps` and `qs` the scattered field at the collocation points of the field mesh is found as

	pf = Mf * ps - Lf * qs;

The complex valued, \f$ F \times 1 \f$ type vector contains the scattered field at the collocation points of the field mesh.

Example {#app_helmholtz_3d_coll_mex_example}
=======

