NiHu 1.1 {#mainpage}
========

[whatsnew11]: \ref news_whatsnew_11 "What's new in release 1.1"

[TOC]

News {#main_news}
====

- **NiHu 1.1 is released** on 2014.03.28.
	- **Support for vector PDE**
	- **Support for the easy incorporation of fundamental solutions**
	- **A generic approach for the accurate evaluation of strongly and hypersingular integrals**
	- For a detailed description of new features visit the \subpage news_whatsnew_11 page

- **Nightly builds** are available from 2013.09.27 built daily at 1 a.m. CET.

- **NiHu 1.0 is released** on 2013.09.13. featuring
	- **C++ template library**
		- A general core functionality for the evaluation of weighted residuals
		- Specialisations for Helmholtz and Laplace kernels
	- **Matlab interface** providing
		- Parametric mesh generation
		- Easy plot generation and post-processing
		- Compatibility with other software
	- **Various tutorials**
		- Guiding you through the implementation of different problems from basic to more complicated ones using the C++ template library of NiHu.
		- Discussing the usage of the Matlab interface for solving different types of boundary element problems by means of implementing mex functions and calling them from Matlab. 
		- Helping you starting the development of your own codes relying on NiHu.

Contact information {#main_contact}
===================

NiHu is developed by the coworkers of the [Laboratory of Acoustics and Studio Technologies](http://last.hit.bme.hu/) at the [Department of Networked Systems and Services](http://www.hit.bme.hu/) of the [Budapest University of Technology and Economics](http://www.bme.hu/).

Members of the development team are:

- [Péter Fiala](mailto:fiala@hit.bme.hu)
- [Péter Rucz](mailto:rucz@hit.bme.hu)

Please feel free to contact any of the developers with your questions or comments regarding NiHu. Your feedbacks are also greatly appreciated.

Licensing {#main_licensing}
=========

NiHu is distributed under the [GNU General Public License](http://www.gnu.org/licenses/gpl.html)

\defgroup core The software core

\defgroup funcspace Function space representations
\ingroup core

\defgroup quadrature Numerical integration
\ingroup core

\defgroup kernel Kernel evaluations
\ingroup core

\defgroup assembly Assembling the system matrices
\ingroup core

\defgroup intop Integral operator and Weighted residual representations
\ingroup core

\defgroup tmp Template metaprogramming libraries
\ingroup core

\defgroup interface Interface with other software

\defgroup matlab Matlab interface
\ingroup interface

\defgroup library Libraries for specific fields of application

