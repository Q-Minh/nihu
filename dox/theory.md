Theory {#theory}
======

[TOC]

Weighted Residuals in a Nuthshell {#weightedresidual}
=================================

From a computational point of view, BEM can be considered as the evaluation of a _Weighted Residual_ \f$R\f$ of the form
\f[ R = \int \int M(x) K(x, y) N(y) dy dx \f]
where \f$K(x,y)\f$ denotes a _Kernel_, \f$M(x)\f$ is a _test function_ and \f$N(y)\f$  is a _trial function_.
The test and trial functions are elements of the _test function space_ and _trial function space_. The function spaces are described by their basis functions \f$M_i(x)\f$ and \f$N_j(y)\f$:
\f[ M(x) = \sum_{i} m_i M_i(x), \quad N(y) = \sum_{j} n_j N_j(x) \f]
So evaluation of the weighed residual is derived from to the evaluation of weighted residuals on the basis functions:
\f[ R_{ij} = \int \int M_i(x) K(x, y) N_j(y) dy dx \f]


Collocational and Galerkin BEM {#collocgalerkin}
------------------------------

When following the _collocational BEM_ approach, the test function is a superposition of Dirac delta functions located at \f$x_i\f$, resulting in the weighted residual expression
\f[ M_i(x) = \delta(x-x_i) \rightarrow R_{ij} = \int K(x_i, y) N_j(y) dy \f]

When following the _Galerkin BEM_ approach, the trial and test functions spaces are identical:
\f[ M_i(x) = N_i(x) \rightarrow R_{ij} = \int \int N_i(x) K(x, y) N_j(y) dy dx \f]


Elements and Fields {#elemfield}
-------------------

In a typical implementation with a general domain geometry, the function spaces' base functions are piecewise polynomial functions over sufficiently small subdomains of the problem geometry:
\f[ R_{ij} = \int_{E_i} \int_{E_j} M_i(x) K(x, y) N_j(y) dy dx \f]
The subdomain is termed _element_. An element is described by its corner nodes \f$x_k\f$ and geometrical interpolation functions \f$L_k(\xi)\f$ as
\f[ E_i = \left\{x : x = \sum_{k} L_k(\xi) x_k, \xi \in \mathcal{D}\right\} \f]
where \f$\mathcal{D}\f$ denotes the element's _base domain_.

An element describes only small scale geometrical behaviour of the integration domain. The interpolation function provides a method to express any internal point \f$x\f$ of an element, or to express the normal vector at any internal point of the element domain.

A set of elements is the problem geometry itself, called _Mesh_.

A single element with the accompaining base functions (\f$M_i(x)\f$ or \f$N_j(x)\f$) is termed a _field_. The set of all fields (that is a mesh with its test or trial basis functions) is the _function space_ itself.

A field is called isoparametric if its trial or test functions are identical with its element's geometrical interpolation functions. A field is called constant if its trial or test functions are constant.

