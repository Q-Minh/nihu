Customising singular integrals {#tut_custom_singular_integrals}
==============================

\page tut_custom_singular_integrals

[BEM example]:\ref theo_bem_example
[operator notation]:\ref bem_example_op
[weighted residuals]:\ref bem_example_res
[Rayleigh integral tutorial]:\ref tut_rayleigh_integral
[function space views]:\ref tut_funcspace_const_iso_view
[dirac view]:\ref tut_funcspace_dirac_view
[TOC]

Introduction {#tut_custom_singular_integrals_intro}
============

This tutorial explains how to customise the evaluation of a singular integral in NiHu. Our examples of demonstration are the collocational singular integrals of the 3D Helmholtz kernels on constant triangle elements.

Theory {#tut_custom_singular_integrals_theory}
======

The singular integrals {#tut_custom_singular_integrals_integrals}
----------------------

The collocational singular integral of the single layer potential kernel on a constant triangular element reads as

\f$
\displaystyle
L = \int_S \frac{\exp(-i k r)}{4\pi r}\mathrm{d}S_{\bf y}, \quad r = |{\bf y} - {\bf x}_0|
\f$

where \f$ {\bf x}_0 \f$ is the singular collocation point in the center of the element.

The singular integral of the double layer potential kernel on a plane element is zero, as the distance vector from the singular point \f$ {\bf x}_0 \f$ to the integral location \f$ {\bf y} \f$ is perpendicular to the element normal \f$ {\bf n}_{\bf y} \f$.

\f$ \displaystyle
M = \int_S \frac{\partial}{\partial n_{\bf y}} \frac{\exp(-i k r)}{4\pi r}\mathrm{d}S_{\bf y} = 0
\f$

For similar reasons, the integral of the transpose of the double layer potential kernel is zero too.

The integral of the hypersingular kernel on a plane triangle element reads as

\f$ \displaystyle
N = \frac{\partial}{\partial n_{\bf x}} \int_S \frac{\partial}{\partial n_{\bf y}} \frac{\exp(-i k r)}{4\pi r}\mathrm{d}S_{\bf y}
\f$

The method of static subtraction {#tut_custom_singular_integrals_subtraction}
--------------------------------

The integral \f$ L \f$ can be regularised by subtracting and adding its static \f$ k = 0 \f$ part:

\f$
\displaystyle
L =
L_0 + L_k =
\int_S \frac{1}{4\pi r}\mathrm{d}S_{\bf y} +
\int_S \frac{\exp(-ikr)-1}{4\pi r}\mathrm{d}S_{\bf y} =
\int_S \frac{1}{4\pi r}\mathrm{d}S_{\bf y} -
\frac{ik}{4\pi} \int_S \exp(-ikr/2) \mathrm{sinc}({kr/2}) \mathrm{d}S_{\bf y}
\f$

The dynamic part is clearly regular and can be integrated with a low order regular quadrature. The static singular part can be evaluated analytically:

\f$
\displaystyle
L_0 = \frac{1}{4\pi}\sum_{i=1}^{3} R_i \sin \alpha_i \left[ \log\left(\frac{\tan{\frac{\phi_1+\alpha_i}{2}}}{\tan{\frac{\alpha_i}{2}}}\right) \right]
\f$


