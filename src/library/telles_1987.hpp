/**
 * @file telles_1987.hpp 
 * @brief Implementation of Telles' quadrature transform method
 * @ingroup lib_singular
 */

#ifndef NIHU_TELLES_1987_HPP_INCLUDED
#define NIHU_TELLES_1987_HPP_INCLUDED

#include "core/quadrature.hpp"

#include <cmath>

namespace NiHu 
{

/**
 * @brief Telles third order polynomial transform 
 * @tparam QuadDerived Derived quadrature type 
 * @param[in] q Quadrature to transform
 * @param[in] eta_bar The singular location in the element coordinate system
 * @param[in] d Distance of the singular point in the other direction
 * @return Transformed quadrature with new base points and weights
 * 
 * @details
 * Performs the third order polynomial transform on a quadrature to cancel
 * singularities or calculate nearly singular integrals more accurately.
 * The method is implemented as given in Ref. [1].
 * 
 * \par
 * <b>References:</b><br>
 *	[1] J.C.F. Telles: <i>A self-adaptive co-ordinate transformation for
 *	    efficient numerical evaluation of general boundary element
 *	    integrals.</i> International Journal for numerical methods in
 *	    engineering. Vol. <b>24</b>, pp. 959-973 (1987)
 *	    DOI: <a href="http://dx.doi.org/10.1002/nme.1620240509">10.1002/nme.1620240509</a>
 */
template <class QuadDerived>
QuadDerived telles_transform(
	quadrature_base<QuadDerived> const &q, 
	typename quadrature_traits<QuadDerived>::domain_t::xi_t const &eta_bar,
	typename quadrature_traits<QuadDerived>::domain_t::scalar_t const &d = typename quadrature_traits<QuadDerived>::domain_t::scalar_t())
{
	typedef quadrature_traits<QuadDerived> traits_t;
	typedef typename traits_t::domain_t domain_t;
	int const n_dims = domain_t::dimension;
	typedef typename domain_t::xi_t xi_t;
	typedef typename domain_t::scalar_t scalar_t;
	
	QuadDerived result = q.derived();
	
	/* Distance dependent parametrization (see (27) in [1] */
	scalar_t r_bar;
	if (d < 0.05)
		r_bar = 0;
	else if (d <= 1.3)
		r_bar = 0.85 + 0.24 * std::log(d);
	else if (d <= 3.618)
		r_bar = 0.893 + 0.0832 * std::log(d);
	else
		r_bar = 1;
	
	/* Go through the dimensions */
	for (int i_dim = 0; i_dim < n_dims; ++i_dim) 
	{
		scalar_t const& e = eta_bar(i_dim);         // Shortcut notation 
		scalar_t f = (1 + 2*r_bar);          		// Helper quantity
		// Intermediate quantities defined in (21)
		scalar_t q = 1/(2*f)*((e*(3-2*r_bar) - 2*e*e*e/f)*1/f - e);
		scalar_t p = 1/(3*f*f)*(4*r_bar*(1-r_bar) + 3*(1 - e*e));
		scalar_t s = std::sqrt(q*q + p*p*p);
		
		// The singular point in gamma domain
		scalar_t gamma_bar = std::cbrt(-q+s) + std::cbrt(-q-s) + e/f;
    
		scalar_t Q = 1 + 3*gamma_bar*gamma_bar;
		// Polynomial coefficients as defined in eq. (21)
		scalar_t a = (1 - r_bar) / Q;
		scalar_t b = -3*(1 - r_bar) * gamma_bar / Q;
		scalar_t c = (3*gamma_bar*gamma_bar + r_bar) / Q;
		scalar_t d = -b;

		/* Go through all quadrature points and scale them */
		for (auto it = result.begin(); it != result.end(); ++it) 
		{
			xi_t &xi = it->get_xi();
			scalar_t x = xi(i_dim);
			/* Update the base point */
			xi(i_dim) = a*x*x*x + b*x*x + c*x + d;
			/* Evaluate jacobian */
			scalar_t jac = 3*a*x*x + 2*b*x + c; 
			/* Update weight */
			it->set_w(it->get_w() * jac);
		}
	}
	
	
	return result;
}

} // end of namespace NiHu

#endif /* NIHU_TELLES_1987_HPP_INCLUDED */
