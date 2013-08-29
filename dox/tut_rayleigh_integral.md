Rayleigh integral {#tut_rayleigh_integral}
=================

\page tut_rayleigh_integral

[TOC]

Introduction {#tut_rayleigh_intro}
============

This tutorial explains how to use NiHu to evaluate the Rayleigh integral.
In this tutorial, we learn how to connect NiHu with Matlab.
The problem geometry and parameters will be imported from, and
the results will be passed back to Matlab using NiHu's Matlab interface.

The Rayleigh integral
---------------------

The Rayleigh integral is used to compute the acoustic pressure \f$ p({\bf x}) \f$ radiated from a vibrating planar radiator surface \f$ S \f$ embedded into an infinite rigid plane.
The vibrating radiator surface is characterised by its normal velocity function \f$ v({\bf y}) \f$.
The radiated pressure is evaluated as

\f$ p({\bf x}) = \int_S G({\bf x}, {\bf y}) v({\bf y}) dS_y \f$

where \f$ G \f$ denotes the fundamental solution of the Helmholtz equation in 3D:

\f$ G({\bf x}, {\bf y}) = \exp(-ikr)/4\pi r, \quad r = |{\bf y} - {\bf x}|\f$


The C++ code {#tut_rayleigh_Cpp_code}
============

Included modules and typedefs {#tut_rayleigh_include}
-----------------------------

\snippet rayleigh_integral_3d.mex.cpp Header

The MEX function
----------------

\snippet rayleigh_integral_3d.mex.cpp Mex function

Mesh generation {#tut_rayleigh_mesh_generation}
---------------

\snippet rayleigh_integral_3d.mex.cpp Surface mesh

\snippet rayleigh_integral_3d.mex.cpp Field point mesh

Function space generation {#tut_rayleigh_function_space}
-------------------------

\snippet rayleigh_integral_3d.mex.cpp Kernel and weighted residual

Kernel and weighted residual {#tut_rayleigh_kernel}
----------------------------


\snippet rayleigh_integral_3d.mex.cpp Kernel and weighted residual


The complete source of the tutorial is found here: tutorial/rayleigh_integral_3d.mex.cpp

