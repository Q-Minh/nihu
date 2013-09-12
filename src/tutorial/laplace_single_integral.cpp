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

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
//! [Mesh]
	dMatrix nodes(4,3);
	nodes << 0., 0., 0., /*|*/ 1., 0., 0., /*|*/ 1., 1., 0., /*|*/ 0., 1., 0.;
	uMatrix elements(1, 1+4);
	elements <<	quad_1_elem::id, 0, 1, 2, 3;
	auto msh = create_mesh(nodes, elements, _quad_1_tag());
//! [Mesh]

//! [Function spaces]
	auto const &trial = constant_view(msh);
	auto const &test = dirac(trial);
//! [Function spaces]

//! [Weighted residual]
	dMatrix A(1, 1);
	A.setZero();

	auto K = create_integral_operator(laplace_3d_SLP_kernel());
	A << ( test * K[trial] );
//! [Weighted residual]

//! [Results]
	std::cout << "WR matrix: " << A << std::endl;

	double anal = std::log(1.+std::sqrt(2.)) / M_PI;
	std::cout << "log10 error = " << std::log10(std::abs(A.sum() / anal - 1.)) << std::endl;
//! [Results]

	return 0;
}

