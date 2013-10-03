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

#include "util/eigen_utils.hpp"
#include "core/weighted_residual.hpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;

template <class test_space_t, class trial_space_t>
void test(test_space_t const &test_space, trial_space_t const &trial_space)
{
	auto I = identity_integral_operator();

	unsigned m = test_space.get_num_dofs();
	unsigned n = trial_space.get_num_dofs();

	dMatrix result(m, n);
	result.setZero();

	result << ( test_space * I[trial_space] );

	std::cout << result << std::endl;
	std::cout << "matrix sum: " << result.sum() << std::endl;
}


template <class mesh_t>
void isoparam_galerkin_test(mesh_t const &mesh)
{
	auto const &f_sp = isoparametric_view(mesh);
	test(f_sp, f_sp);
}

template <class mesh_t>
void constant_galerkin_test(mesh_t const &mesh)
{
	auto const &f_sp = constant_view(mesh);
	test(f_sp, f_sp);
}

template <class mesh_t>
void constant_collocational_test(mesh_t const &mesh)
{
	auto const &f_sp = constant_view(mesh);
	test(dirac(f_sp), f_sp);
}


int main(void)
{
	dMatrix nodes(9,3);
	uMatrix elements(4,5);

	double cx = 2.0, cy = 3.0;

    nodes <<
		-cx, -cy, 0.,
		0., -cy, 0.,
		cx, -cy, 0.,

		-cx, 0., 0.,
		0., 0., 0.,
		cx, 0., 0.,

		-cx, cy, 0.,
		0., cy, 0.,
		cx, cy, 0.;

    elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		quad_1_elem::id, 1, 2, 5, 4,
		quad_1_elem::id, 3, 4, 7, 6,
		quad_1_elem::id, 4, 5, 8, 7;


	auto mesh = create_mesh(nodes, elements, quad_1_tag());

	std::cout << std::endl << "testing isoparametric galerkin case" << std::endl;
	isoparam_galerkin_test(mesh);

	std::cout << std::endl << "testing constant galerkin case" << std::endl;
	constant_galerkin_test(mesh);

	std::cout << std::endl << "testing constant collocational case" << std::endl;
	constant_collocational_test(mesh);

	return 0;
}
