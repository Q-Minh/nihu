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
#include "library/quad_1_gauss_field.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

int main(void)
{
	std::cout << NiHu::quad_1_gauss_field::elem_t::domain_t::dimension << std::endl;
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

	uMatrix fields(4, 1+4+4);
	fields <<
		NiHu::quad_1_gauss_field::id, 0, 1, 4, 3,  0, 1, 2, 3,
		NiHu::quad_1_gauss_field::id, 1, 2, 5, 4,  4, 5, 6, 7,
		NiHu::quad_1_gauss_field::id, 3, 4, 7, 6,  8, 9, 10, 11,
		NiHu::quad_1_gauss_field::id, 4, 5, 8, 7,  12, 13, 14, 15;

	// create function space using the factory function
	auto fsp = NiHu::create_function_space(nodes, fields, NiHu::quad_1_gauss_field_tag());

	// allocate and clear the result matrix
	dMatrix A(fsp.get_num_dofs(), fsp.get_num_dofs());
	A.setZero();

	// and evaluate the weighed residual into our result matrix
	auto L = NiHu::create_integral_operator(NiHu::laplace_3d_SLP_kernel());
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


