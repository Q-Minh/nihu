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
L_0 = \frac{1}{4\pi}\sum_{i=1}^{3} R_i \sin \alpha_i \left[ \log\left(\frac{\tan{\frac{\phi_i+\alpha_i}{2}}}{\tan{\frac{\alpha_i}{2}}}\right) \right]
\f$

where $R_i$ denotes the distance from the singular point to the \f$ i \f$-th corner of the triangle, and the angles \f$ \phi_i \f$ and \f$ \alpha_i \f$ are explained in the figure below:

The C++ code {#tut_custom_singular_integrals_cpp}
============

First we define a helper function `triangle_helper` that computes the radii \f$ R_i \f$ and the angles \f$ \alpha_i \f$ and \f$ \theta_i \f$ for a triangle element. This function will be used for the evaluation of both singular integrals:
~~~~~~~~~~~
void triangle_helper(
	tria_1_elem const &elem,
	double r[],			// output radii
	double theta[],		// output angles
	double alpha [])	// output angles
{
	auto const &C_old = elem.get_coords();		// element vertices in columns
	auto const &x0 = elem.get_center();			// collocational point x0

	typename tria_1_elem::coords_t R, C;
	for (unsigned i = 0; i < 3; ++i)
	{
		R.col(i) = C_old.col(i) - x0;	// vector from x0 to corner
		r[i] = R.col(i).norm();			// distance
		R.col(i) /= r[i];				// normalise vector
		C.col(i) = C_old.col(i) - C_old.col((i+1) % 3);	// side vector
		C.col(i) /= C.col(i).norm();					// normalised
	}

	for (unsigned i = 0; i < 3; ++i)
	{
		theta[i] = std::acos(R.col(i).dot(R.col((i+1) % 3)));	// angle output
		alpha[i] = std::acos(R.col(i).dot(C.col(i)));			// angle output
	}
}
~~~~~~~~~~~

The singular integral of the single layer potential kernel for the constant triangle element is customised by specialising the class template ::singular_integral_shortcut for the specified singular integral type.

The class template is declared as
~~~~~~~~~~~~
template <class Kernel, class TestField, class TrialField, class Singulartiy, class Enable>
class singular_integral_shortcut;
~~~~~~~~~~~~
Where the parameters are
- the kernel type (::helmholtz_3d_SLP_kernel in our case)
- the test and trial field types (both constant triangles, the test field is a Dirac-view)
- the singularity type (::singularity::face_match_type in our case)
- an additional parameter to make complex type selection easy using std::enable_if


