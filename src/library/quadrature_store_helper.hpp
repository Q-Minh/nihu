#ifndef QUADRATURE_STORE_HELPER_HPP_INCLUDED
#define QUADRATURE_STORE_HELPER_HPP_INCLUDED

#include "../core/quadrature.hpp"


namespace NiHu
{

/** \brief store-wrapper of a statically stored quadrature */
template <class domain_t, size_t order>
struct regular_quad_store
{
	/** \brief the stored static quadrature member */
	static gaussian_quadrature<domain_t> const quadrature;
};

/** \brief definition of the statically stored quadrature member */
template <class domain_t, size_t order>
gaussian_quadrature<domain_t> const regular_quad_store<domain_t, order>::quadrature(order);

}

#endif
