// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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

//! [Forward declaration]
class quad_1_gauss_shape_set;
//! [Forward declaration]

//! [Shape traits]
template<>
struct shape_set_traits<quad_1_gauss_shape_set>
{
	typedef quad_domain domain_t;
	static unsigned const num_nodes = 4;
	static unsigned const polynomial_order = 1;
	static unsigned const jacobian_order = 1;
};
//! [Shape traits]


//! [Shape class]
class quad_1_gauss_shape_set
	: public shape_set_base<quad_1_gauss_shape_set>
{
public:
	static shape_t eval_shape(xi_t const &_xi);
	static dshape_t eval_dshape(xi_t const & _xi);
	static xi_t const *corner_begin(void);

protected:
	static const xi_t m_corners[num_nodes];
};
//! [Shape class]

//! [Shape lsets]
inline quad_1_gauss_shape_set::shape_t
	quad_1_gauss_shape_set::eval_shape(quad_1_gauss_shape_set::xi_t const &_xi)
{
	scalar_t xi = _xi[0], eta = _xi[1];
	shape_t L;
	L <<
		(1.0 - std::sqrt(3.0)*xi) * (1.0 - std::sqrt(3.0)*eta) / 4.0,
		(1.0 + std::sqrt(3.0)*xi) * (1.0 - std::sqrt(3.0)*eta) / 4.0,
		(1.0 - std::sqrt(3.0)*xi) * (1.0 + std::sqrt(3.0)*eta) / 4.0,
		(1.0 + std::sqrt(3.0)*xi) * (1.0 + std::sqrt(3.0)*eta) / 4.0;
	return L;
}

inline quad_1_gauss_shape_set::dshape_t
	quad_1_gauss_shape_set::eval_dshape(quad_1_gauss_shape_set::xi_t const & _xi)
{

	scalar_t xi = _xi[0], eta = _xi[1];

	dshape_t dL;
	dL <<
		-std::sqrt(3.0) * (1.0 - std::sqrt(3.0)*eta) / 4.0, (1.0 - std::sqrt(3.0)*xi) * -std::sqrt(3.0) / 4.0,
		+std::sqrt(3.0) * (1.0 - std::sqrt(3.0)*eta) / 4.0, (1.0 + std::sqrt(3.0)*xi) * -std::sqrt(3.0) / 4.0,
		-std::sqrt(3.0) * (1.0 + std::sqrt(3.0)*eta) / 4.0, (1.0 - std::sqrt(3.0)*xi) * +std::sqrt(3.0) / 4.0,
		+std::sqrt(3.0) * (1.0 + std::sqrt(3.0)*eta) / 4.0, (1.0 + std::sqrt(3.0)*xi) * +std::sqrt(3.0) / 4.0;
	return dL;
}
//! [Shape lsets]

//! [Shape corners]
inline quad_1_gauss_shape_set::xi_t const *
	quad_1_gauss_shape_set::corner_begin(void)
{
	return m_corners;
}

quad_1_gauss_shape_set::xi_t
	const quad_1_gauss_shape_set::m_corners[quad_1_gauss_shape_set::num_nodes] = {
		quad_1_gauss_shape_set::xi_t(-std::sqrt(3.0)/3.0, -std::sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(+std::sqrt(3.0)/3.0, -std::sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(-std::sqrt(3.0)/3.0, +std::sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(+std::sqrt(3.0)/3.0, +std::sqrt(3.0)/3.0)
};
//! [Shape corners]


//! [Field typedef]
typedef field<quad_1_elem, quad_1_gauss_shape_set> quad_1_gauss_field;
//! [Field typedef]

//! [Field id]
template <>
struct field_id<quad_1_gauss_field>
{
	static unsigned const value = 666;
};
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


