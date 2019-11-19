Plane wave scattering from a rigid sphere {#bench_plane_wave_sphere}
=========================================

\page bench_plane_wave_sphere Plane wave scattering from a rigid sphere

[TOC]

[ref1]:http://ansol.us/Products/Coustyx/Validation/MultiDomain/Scattering/PlaneWave/HardSphere/Downloads/dataset_description.pdf

Problem definition {#bench_plane_wave_sphere_def}
==================

The incident field is given as

\f$ 
    \displaystyle p_{\mathrm{inc}} = \mathrm{e}^{-\mathrm{j} \mathbf{k} \cdot \mathbf{x}}
\f$

The incident field can be expanded into a series of spherical functions in the \f$ (r,\theta) \f$ coordinate frame as

\f$
    \displaystyle p_{\mathrm{inc}}(r,\theta) = \sum_{l=1}^{\infty} (2l+1) \mathrm{j}^{l} j_{l}^{(1)} (kr) P_l(\cos(\theta))
\f$

The scatterer is perfectly rigid, i.e., the normal derivative of the sound pressure is zero on its surface
\f$
    \displaystyle q = p'_n = 0
\f$

\f$
    \displaystyle q_{\mathrm{scat}} = -q_{\mathrm{inc}}
\f$

The scattered field is found by satisfying the boundary condition and using the spherical expansion of the incident field

\f$
    \displaystyle p_{\mathrm{scat}}(r, \theta) = \sum_{l=0}^{\infty} A_l h_l^{(2)}(kr) P_l(\cos(\theta))
\f$

The coefficient \f$ A_l \f$ is found as

\f$
    \displaystyle A_l = -(2l + 1) \mathrm{j}^l \frac{l j_{l-1}^{(1)}(kr_0) - (l+1)j_{l+1}^{(1)}(kr_0)}{lh_{l-1}^{(2)}(kr_0) - (l+1)h_{l+1}^{(2)}(kr_0)}
\f$

Finally, the total sound pressure field is attained as the sum of the incident and scattered fields as

\f$
    \displaystyle p(r, \theta) = p_{\mathrm{inc}}(r, \theta) + p_{\mathrm{scat}}(r, \theta)
\f$

BEM solution {#bench_plane_wave_sphere_sol}
============

The BEM solution of the above problem is found by prescribing the scattered velocity field

The infinite sums in series expansions of the incident and scattered fields are truncated up to a given order.
To ensure that the error introduced by the truncation is small, the truncated series expansion and the plane wave formula for the incident field are compared and the relative error is evaluated.
It is found that it is enough to keep a few tens of terms in the expansion to get a relative error of \f$ \varepsilon \approx 10^{-7} \f$.

Results {#bench_planw_wave_sphere_res}
=======

    Series truncation error:  1.29315e-07
    Surface relative error:   0.0140948
    Field relative error:     0.013953

The figure below shows the absolute 
    
\image html bench_plane_wave_sphere.png
