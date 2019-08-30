Fast Multipole Method for the 2D Helmholtz Kernel {#tut_fmm_helmholtz2d}
=================================================

\page tut_fmm_helmholtz2d Evaluating the Helmholtz kernel using the Fast Multipole Method

[TOC]

Introduction {#tut_fmmh2d_intro}
============

In this example we compute the acoustic field radiated by a set of 2D acoustic 
monopoles.
This problem involves the evaluation of the 2d Helmholtz kernel between source and
receiver locations and adding the sources' contributions.
We solve this matrix multiplication by means of the wideband fast multipole method
for the 2D Helmholtz kernel.


Code {#tut_fmmh2d_code}
====

The below sections present the full code needed to compute the MVP (matrix vector product).

Type Definitions {#tut_fmmh2d_types}
----------------

A few types are defined first fo convenience, including the cluster type,
 cluster tree type and physical location type.

\snippet fmm_helmholtz2d.cpp typedefs

Source and receiver nodes  {#tut_fmmh2d_nodes}
-------------------------

We define random source and receiver locations:

\snippet fmm_helmholtz2d.cpp source and receiver nodes

Cluster tree  {#tut_fmmh2d_tree}
------------

A cluster tree is built that contains all the source and receiver nodes.
The clusters are divided until they contain less than 10 nodes.

\snippet fmm_helmholtz2d.cpp tree and lists

FMM initializaion {#tut_fmmh2d_fmm}
-----------------

We initialize further data structures needed to compute the fmm operators.
This is done by initializing the fmm object, and initializing its level data vector using the cluster tree.
The level data vector is attached to the tree by linking the appropriate level data structures to each cluster.


\snippet fmm_helmholtz2d.cpp fmm initialization

FMM operator manipulations  {#tut_fmmh2d_operators}
--------------------------

We get the fmm operators from the fmm object, and combine them into a collection for later convenience.
Two functors are defined, first to assign the nodes and the cluster tree to the operators, and a second to assign the interaction lists to them.

Using the two functors, all the fmm operators are easily indexed and precomputed.


\snippet fmm_helmholtz2d.cpp operator manipulations

FMM matrix and MVP  {#tut_fmmh2d_matrix}
------------------

Finally, an fmm matrix is instantiated from the precomputed operators, a source vector is allocated, and the response is computed by simple matrix vector product.

\snippet fmm_helmholtz2d.cpp matrix
