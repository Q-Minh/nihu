#ifndef LINE_QUAD_STORE_HPP_INCLUDED
#define LINE_QUAD_STORE_HPP_INCLUDED

#include "../core/gaussian_quadrature.hpp"

namespace NiHu
{

/** \brief store-wrapper of a statically stored line quadrature */
template <unsigned order>
struct line_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<line_domain> const quadrature;
};

/** \brief definition of the statically stored line quadrature member */
template <unsigned order>
gaussian_quadrature<line_domain> const line_quad_store<order>::quadrature(order);

} // end of namespace NiHu

#endif

