/** \file line_1_gauss_shape_set.hpp 
 * \brief implementation of shape set \ref NiHu::line_1_gauss_shape_set
 */

#ifndef LINE_1_GAUSS_SHAPE_SET_HPP_INCLUDED
#define LINE_1_GAUSS_SHAPE_SET_HPP_INCLUDED

#include "../core/shapeset.hpp"
#include "lib_domain.hpp"

// define the new shape set's traits
namespace NiHu
{
	
// predefine the new shape set
class line_1_gauss_shape_set;

namespace shape_set_traits
{
	/** \todo name is not defined */

	// the shape set is defined over the line domain
	template <>
	struct domain<line_1_gauss_shape_set> : line_domain {};

	// it has 2 shapeset nodes
	template <>
	struct num_nodes<line_1_gauss_shape_set> { enum { value = 2 }; };

	/** \todo id is not defined */

	// N(xi) is first order
	template <>
	struct polynomial_order<line_1_gauss_shape_set> { enum { value = 1 }; };

	// if used as geometrical shape set, Jacobian is 0-st order in xi
	template <>
	struct jacobian_order<line_1_gauss_shape_set> { enum { value = 0 }; };

	// N(xi) is computed on the fly
	template <>
	struct shape_complexity<line_1_gauss_shape_set, 0>
	{
		typedef matrix_function_complexity::general type;
	};

	// N'(xi) is constant
	template <>
	struct shape_complexity<line_1_gauss_shape_set, 1>
	{
		typedef matrix_function_complexity::constant type;
	};

	// N''(xi) is zero
	template <>
	struct shape_complexity<line_1_gauss_shape_set, 2>
	{
		typedef matrix_function_complexity::zero type;
	};

	// DOF locations are inside the element (1DOF)
	template <>
	struct position_dof_vector<line_1_gauss_shape_set>
	{
		typedef tmp::vector<dof1, dof1> type;
	};
} // end of namespace shape_set_traits


// define the new shape set class
class line_1_gauss_shape_set
	: public NiHu::shape_set_base<line_1_gauss_shape_set>
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

// define expression for computing the zeroth order derivative N(xi)
template <>
class shape_function<line_1_gauss_shape_set, 0>
{
	typedef shape_set_traits::shape_value_type<line_1_gauss_shape_set, 0>::type shape_t;
	typedef shape_set_traits::domain<line_1_gauss_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &_xi)
	{
		auto xi = _xi[0];
		return ( line_1_gauss_shape_set::shape_t() <<
			(1.0 - std::sqrt(3.0)*xi),
			(1.0 + std::sqrt(3.0)*xi)
		).finished() / 2.0;
	}
};


// define expression for computing the first order derivative N(xi)
template <>
class shape_function<line_1_gauss_shape_set, 1>
{
	typedef shape_set_traits::shape_value_type<line_1_gauss_shape_set, 1>::type shape_t;
	typedef shape_set_traits::domain<line_1_gauss_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return ( line_1_gauss_shape_set::dshape_t() <<
			-std::sqrt(3.0),
			+std::sqrt(3.0)
		).finished() / 2.0;
	}
};


// define expression for computing the first order derivative N(xi)
template <>
class shape_function<line_1_gauss_shape_set, 2>
{
	typedef shape_set_traits::shape_value_type<line_1_gauss_shape_set, 2>::type shape_t;
	typedef shape_set_traits::domain<line_1_gauss_shape_set>::type::xi_t xi_t;
public:
	static shape_t eval(xi_t const &)
	{
		return shape_t::Zero();
	}
};

} // end of namespace NiHu

#endif // LINE_1_GAUSS_SHAPE_SET_HPP_INCLUDED
