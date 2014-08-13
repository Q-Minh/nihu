// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
#include "library/lib_element.hpp"

//! [Forward declaration]
class quad_1_gauss_shape_set;
//! [Forward declaration]

//! [Shape traits]
namespace shape_set_traits
{
	template <>
	struct domain<quad_1_gauss_shape_set> : quad_domain {};

	template <>
	struct num_nodes<quad_1_gauss_shape_set>
	{
		enum { value = 4 };
	};

	template <>
	struct jacobian_order<quad_1_gauss_shape_set>
	{
		enum { value = 1 };
	};

	template <>
	struct polynomial_order<quad_1_gauss_shape_set>
	{
		enum { value = 1 };
	};

	template <unsigned Order>
	struct shape_complexity<quad_1_gauss_shape_set, Order>
	{
		typedef matrix_function_complexity::general type;
	};

	template <>
	struct position_dof_vector<quad_1_gauss_shape_set>
	{
		typedef tmp::vector<dof2, dof2, dof2, dof2> type;
	};
}
//! [Shape traits]


//! [Shape class]
class quad_1_gauss_shape_set
	: public shape_set_base<quad_1_gauss_shape_set>
{
public:
	static xi_t const *corner_begin(void)
	{
		return m_corners;
	}

protected:
	static const xi_t m_corners[num_nodes];
};

quad_1_gauss_shape_set::xi_t
	const quad_1_gauss_shape_set::m_corners[quad_1_gauss_shape_set::num_nodes] = {
		quad_1_gauss_shape_set::xi_t(-std::sqrt(3.0)/3.0, -std::sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(+std::sqrt(3.0)/3.0, -std::sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(-std::sqrt(3.0)/3.0, +std::sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(+std::sqrt(3.0)/3.0, +std::sqrt(3.0)/3.0)
};
//! [Shape class]

//! [Shape lsets]
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
//! [Shape lsets]

//! [Field typedef]
typedef field<quad_1_elem, quad_1_gauss_shape_set> quad_1_gauss_field;
//! [Field typedef]

//! [Field id]
namespace field_traits
{
	template <>
	struct id<quad_1_gauss_field>
	{
		enum { value = 666 };
	};
}
//! [Field id]

//! [Field tag]
struct gauss_field_tag {};

template <>
struct tag2field<gauss_field_tag> :
	quad_1_gauss_field {};
//! [Field tag]


//! [main]
// double matrix
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
// unsigned matrix
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	std::cout << quad_1_gauss_field::elem_t::domain_t::dimension << std::endl;
	// nodal coordinates of our mesh
	dMatrix nodes(9,3);
	nodes <<
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		2.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		1.0, 1.0, 0.0,
		2.0, 1.0, 0.0,
		0.0, 2.0, 0.0,
		1.0, 2.0, 0.0,
		2.0, 2.0, 0.0;

	// fields (node indices + DOF indices)
	uMatrix fields(4, 1+4+4);
	fields <<
	//  field id                nodes        dofs
		quad_1_gauss_field::id, 0, 1, 4, 3,  0, 1, 2, 3,
		quad_1_gauss_field::id, 1, 2, 5, 4,  4, 5, 6, 7,
		quad_1_gauss_field::id, 3, 4, 7, 6,  8, 9, 10, 11,
		quad_1_gauss_field::id, 4, 5, 8, 7,  12, 13, 14, 15;

	// create function space using the factory function
	auto fsp = create_function_space(nodes, fields, gauss_field_tag());

	// allocate and clear the result matrix
	dMatrix A(fsp.get_num_dofs(), fsp.get_num_dofs());
	A.setZero();

	// and evaluate the weighed residual into our result matrix
	auto L = create_integral_operator(laplace_3d_SLP_kernel());
	A << fsp * L[fsp];

	// print the result matrix
	std::cout << "matrix coefficients\n" << A << "\n\n";
	// Barzini is happy oleoleole leoleole leoleole!
	std::cout << "matrix sum: " << A.sum() << "\n\n";
	double anal = 8.0/M_PI*(std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0);
	std::cout << "analytical: " << anal << "\n\n";
	std::cout << "log10 error: " << std::log10(std::abs(A.sum()/anal - 1.0)) << "\n\n";

	return 0;
}
//! [main]

