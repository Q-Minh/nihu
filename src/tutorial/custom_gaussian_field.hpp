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

class quad_1_gauss_shape_set;

template<>
struct shape_set_traits<quad_1_gauss_shape_set>
{
	typedef quad_domain domain_t;
	static unsigned const num_nodes = 4;
	static unsigned const polynomial_order = 1;
	static unsigned const jacobian_order = 1;
};


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


typedef field<quad_1_elem, quad_1_gauss_shape_set> quad_1_gauss_field;

template <>
struct field_id<quad_1_gauss_field>
{
	static unsigned const value = 666;
};

struct gauss_field_tag {};

template <>
struct tag2field<gauss_field_tag> :
	quad_1_gauss_field {};


