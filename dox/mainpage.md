NiHu 1.0 {#mainpage}
========

[GPL]:http://www.gnu.org/licenses/gpl.html 
 
 News {#main_news}
=================

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

NiHu is developed by the coworkers of the [Laboratory of Acoustics and Studio Technologies](http://last.hit.bme.hu/) at the [Depertment of Networked Systems and Services](http://www.hit.bme.hu/) of the [Budapest University of Technology and Economics](http://www.bme.hu/).

Members of the development team are:

- Péter Fiala [fiala@hit.bme.hu](mailto:fiala@hit.bme.hu)
- Péter Rucz [rucz@hit.bme.hu](mailto:rucz@hit.bme.hu)

Please feel free to contact any of the developers with your questions or comments regarding NiHu. Your feedbacks are also greatly appreciated.

Licensing {#main_licensing}
=========

NiHu is distributed under the GNU General Public License [GPL].

\defgroup core The software core

\defgroup funcspace Function space representations
\ingroup core

\defgroup quadrature Numerical integration
\ingroup core

\defgroup intop Integral operator and Weighted residual representations
\ingroup core

\defgroup tmp Template metaprogramming libraries
\ingroup core

\defgroup interface Interface with other software

\defgroup matlab Matlab interface
\ingroup interface

\defgroup library Libaries for specific fields of application

