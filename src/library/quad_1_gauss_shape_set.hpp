#ifndef QUAD_1_GAUSS_SHAPE_SET_HPP_INCLUDED
#define QUAD_1_GAUSS_SHAPE_SET_HPP_INCLUDED

#include "../core/shapeset.hpp"
#include "lib_domain.hpp"

// define the new shape set's traits
namespace NiHu
{
	
// predefine the new shape set
class quad_1_gauss_shape_set;

namespace shape_set_traits
{
	// the shape set is defined over the quad domain
	template <>
	struct domain<quad_1_gauss_shape_set> : quad_domain {};

	// it has 4 shapeset nodes
	template <>
	struct num_nodes<quad_1_gauss_shape_set> { enum { value = 4 }; };

	// if used as geometrical shape set, Jacobian is 1-st order in xi/eta
	template <>
	struct jacobian_order<quad_1_gauss_shape_set> { enum { value = 1 }; };

	// N(xi,eta) is first order
	template <>
	struct polynomial_order<quad_1_gauss_shape_set> { enum { value = 1 }; };

	// N(xi, eta) and its derivatives are computed on the fly
	template <unsigned Order>
	struct shape_complexity<quad_1_gauss_shape_set, Order>
	{
		typedef matrix_function_complexity::general type;
	};

	// DOF locations are inside the element (2DOF)
	template <>
	struct position_dof_vector<quad_1_gauss_shape_set>
	{
		typedef tmp::vector<dof2, dof2, dof2, dof2> type;
	};
} // end of namespace shape_set_traits


// define the new shape set class
class quad_1_gauss_shape_set
	: public NiHu::shape_set_base<quad_1_gauss_shape_set>
{
public:
	// return iterator to corners
	static xi_t const *corner_begin_impl(void)
	{
		return m_corners;
	}

protected:
	static const xi_t m_corners[num_nodes];
};

// define expression for computing the zeroth order derivative N(xi, eta)
template <>
class shape_function<quad_1_gauss_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<quad_1_gauss_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<quad_1_gauss_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0], eta = _xi[1];
		return ( quad_1_gauss_shape_set::shape_t() <<
			(1.0 - std::sqrt(3.0)*xi) * (1.0 - std::sqrt(3.0)*eta),
			(1.0 + std::sqrt(3.0)*xi) * (1.0 - std::sqrt(3.0)*eta),
			(1.0 - std::sqrt(3.0)*xi) * (1.0 + std::sqrt(3.0)*eta),
			(1.0 + std::sqrt(3.0)*xi) * (1.0 + std::sqrt(3.0)*eta)
		).finished() / 4.0;
	}
};


// define expression for computing the first order derivative N(xi, eta)
template <>
class shape_function<quad_1_gauss_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<quad_1_gauss_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<quad_1_gauss_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0], eta = _xi[1];
		return ( quad_1_gauss_shape_set::dshape_t() <<
			-std::sqrt(3.0) * (1.0 - std::sqrt(3.0)*eta) , (1.0 - std::sqrt(3.0)*xi) * -std::sqrt(3.0),
			+std::sqrt(3.0) * (1.0 - std::sqrt(3.0)*eta) , (1.0 + std::sqrt(3.0)*xi) * -std::sqrt(3.0),
			-std::sqrt(3.0) * (1.0 + std::sqrt(3.0)*eta) , (1.0 - std::sqrt(3.0)*xi) * +std::sqrt(3.0),
			+std::sqrt(3.0) * (1.0 + std::sqrt(3.0)*eta) , (1.0 + std::sqrt(3.0)*xi) * +std::sqrt(3.0)
		).finished() / 4.0;
	}
};

} // end of namespace NiHu

#endif // QUAD_1_GAUSS_SHAPE_SET_HPP_INCLUDED
