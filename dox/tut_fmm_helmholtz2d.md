Fast Multipole Method for the 2D Helmholtz Kernel {#tut_fmm_helmholtz2d}
=================================================

\page tut_fmm_helmholtz2d Evaluating the Helmholtz kernel using the Fast Multipole Method

[TOC]

Introduction {#tut_fmmh2d_intro}
============

In this example we compute the acoustic field radiated by a set of 2D acoustic 
monopoles.
This problem involves the evaluatioj of the 2d Helmholtz kernel between source and
receiver locations and adding the sources' contributions.
We solve this matrixmultiplication by means of the wideband fast multipole method
for the 2D Helmholtz kernel.


Code {#tut_fmmh2d_code}
====

Type Definitions {#tut_fmmh2d_types}
----------------

\snippet fmm_helmholtz2d.cpp typedefs

Source and receiver nodes  {#tut_fmmh2d_nodes}
-------------------------

\snippet fmm_helmholtz2d.cpp source and receiver nodes

Source and receiver nodes  {#tut_fmmh2d_tree}
-------------------------

\snippet fmm_helmholtz2d.cpp tree and lists

FMM initializaion {#tut_fmmh2d_fmm}
-----------------

\snippet fmm_helmholtz2d.cpp fmm initialization

FMM operator manipulations  {#tut_fmmh2d_operators}
--------------------------

\snippet fmm_helmholtz2d.cpp operator manipulations

FMM matrix and MVP  {#tut_fmmh2d_matrix}
------------------

\snippet fmm_helmholtz2d.cpp matrix
