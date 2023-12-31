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

#include <boost/math/constants/constants.hpp>

#include "nihu/core/weighted_residual.hpp"
#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/lib_element.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	using namespace boost::math::double_constants;
	
//! [Mesh]
	dMatrix nodes(4,3);
	nodes << 0., 0., 0., /*|*/ 1., 0., 0., /*|*/ 1., 1., 0., /*|*/ 0., 1., 0.;
	uMatrix elements(1, 1+4);
	elements <<	NiHu::quad_1_elem::id, 0, 1, 2, 3;
	auto msh = NiHu::create_mesh(nodes, elements, NiHu::quad_1_tag());
//! [Mesh]

//! [Function spaces]
	auto const &trial = NiHu::constant_view(msh);
	auto const &test = NiHu::dirac(trial);
//! [Function spaces]

//! [Weighted residual]
	dMatrix A(1, 1);
	A.setZero();

	auto K = NiHu::create_integral_operator(NiHu::laplace_3d_SLP_kernel());
	A << test * K[trial];
//! [Weighted residual]

//! [Results]
	std::cout << "WR matrix: " << A << std::endl;

	double anal = std::log(1. + root_two) / pi;
	std::cout << "log10 error = " << std::log10(std::abs(A.sum() / anal - 1.)) << std::endl;
//! [Results]

	return 0;
}

