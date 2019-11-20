Circular membrane with static load {#bench_circular_membrane}
==================================

\page bench_circular_membrane Circular membrane with static load

[TOC]

Introduction {#bench_circ_memb_intro}
============

This benchmark test evaluates the displacement of a circular membrane under static load.
The membrane is excited by a force concentrated in one point and the boundary of the membrane is assumed to be at rest.
The resulting displacement field on the membrane is evaluted by means of the BEM and is compared to the analytical solution of the same problem.
The analytical solution is attained utilizing the Green's function of the circular membrane based on the method of mirror images.

Analytical solution
===================

The boundary value problem is defined using the two-dimensional Laplace equation as

\f$
	\displaystyle \nabla^2 u(\mathbf{x}) = -\frac{1}{T} g(\mathbf{x}) \qquad \mathbf{x} \in \Omega : \left| \mathbf{x} \right| \le R 
\f$ 

\f$
	\displaystyle u(\mathbf{x}) = 0 \qquad \mathbf{x} \in \Gamma: \left|\mathbf{x}\right| = R
\f$

with \f$ u \f$ denoting the displacement of the membrane, and \f$ g \f$ representing the excitation that has a dimension of force over unit area.
The tension of the membrane \f$ T \f$ is constant and is taken as unity in the following for the sake of convenience.

The free field Green's function of the 2D Laplace equation \f$ G_0(\mathbf{x}, \mathbf{x}_0) \f$ is given as

\f$
	\displaystyle G_0(\mathbf{x}, \mathbf{x}_0) = -\frac{1}{2 \pi} \log \left| \mathbf{x} - \mathbf{x}_0 \right|
\f$

The free field Green's function in itself does not satisfy the prescribed boundary condition, but a modified Green's function can be constructed using the method of mirror images as follows.
For all points in the interior of the membrane \f$ | \mathbf{x}_0 | < R \f$ there exists an image point 

\f$ 
	\displaystyle \mathbf{x}_0^\star = \mathbf{x}_0 \frac{R^2}{| \mathbf{x}_0 |^2}
\f$ 

such that the equality

\f$
	\displaystyle \left|\mathbf{x} - \mathbf{x}_0 \right|^2 = \frac{|\mathbf{x_0}|^2}{R^2} \left| \mathbf{x}  - \mathbf{x}_0^\star \right|^2 
\f$

holds for all points \f$ |\mathbf{x}| = R \f$.

The above property is utilized for formulating the Green's function \f$ G(\mathbf{x}, \mathbf{x}_0) \f$.
By compensating the free field Green's function not being zero at the boundary, the Green's function for the circular membrane reads as

\f$
	\displaystyle G(\mathbf{x}, \mathbf{x_0}) = \frac{1}{2\pi}\left( -\log|\mathbf{x} - \mathbf{x}_0| + \log|\mathbf{x} - \mathbf{x}_0^\star| + \log\frac{|\mathbf{x}_0|}{R} \right)
\f$

In the present case, concentrated force acts on the membrane, i.e. \f$ g(\mathbf{x}) = \delta(\mathbf{x} - \mathbf{x}_{\mathrm{exc}}) \f$ is assumed.
Then, the analytical solution is found by substituting \f$ \mathbf{x}_{\mathrm{exc}} \f$ as \f$ \mathbf{x}_0 \f$ into the Green's function.
Similarly, in case of a distributed load, the resulting displacement field is found by evaluating the convolution integral

\f$
	\displaystyle u(\mathbf{x}) = \int_{\Omega} G(\mathbf{x}, \mathbf{x}_0) g(\mathbf{x}_0) \, \mathrm{d} \mathbf{x}_0
\f$

Boundary element solution
=========================

The incident displacement field is evaluated using the free field Green's function \f$ G_0(\mathbf{x}, \mathbf{x}_0) \f$ 

\f$ 
	\displaystyle u_{\mathrm{inc}}(\mathbf{x}) = G_0(\mathbf{x}, \mathbf{x}_{\mathrm{exc}})
\f$

Zero displacement on the boundary is satisfied by prescribing the scattered field, 

\f$ 
	\displaystyle u_{\mathrm{scat}}(\mathbf{x}) = -u_{\mathrm{inc}}(\mathbf{x}) \qquad \mathbf{x} \in \Gamma 
\f$

Then, the normal derivative of the scattered displacement field is evaluted at the boundary by solving the BEM system of equations

\f$
	\displaystyle \mathbf{L}_{\mathrm{b}} \mathbf{q}_{\mathrm{b}} = \mathbf{M}_{\mathrm{b}} \mathbf{u}_{\mathrm{b}} 
\f$

Then, the scattered displacement field over the whole membrane is attained by evaluating the matrix-vector product

\f$
	\displaystyle u_{\mathrm{scat}} = \mathbf{L}_{\mathrm{m}} \mathbf{q}_{\mathrm{b}} - \mathbf{M}_{\mathrm{m}} \mathbf{u}_{\mathrm{b}}
\f$

Finally, the total displacement results as the sum of the incident and scattered displacement fields

\f$
	\displaystyle u = u_{\mathrm{inc}} + u_{\mathrm{scat}}
\f$ 

Results {#bench_circ_memb_results}
=======



	Relative error of displacement on membrane: 0.0119008

If the boundary is discretized using quadratic elements, the relative error is remarkably reduced even if the number of elements is small.

	Relative error of displacement on membrane: 0.00185664

The resulting displacement is illustrated in the figure below, with the red cross denoting the point where the excitation is applied.

\image html bench_circular_membrane.png