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

//! [Includes]
#include "nihu/core/weighted_residual.hpp"
#include "nihu/library/laplace_kernel.hpp"
#include "nihu/library/lenoir_salles_2012.hpp"
#include "nihu/library/lib_element.hpp"
//! [Includes]

//! [Typedefs]
// dynamically resizeable double matrix
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
// dynamically resizeable unsigned matrix
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
//! [Typedefs]

int main(void)
{
	using namespace boost::math::double_constants;

	// nodal coordinates in 3D, 3 nodes
	dMatrix nodes(9,3);
	nodes <<
		-1.0, -1.0, 0.0,
		 0.0, -1.0, 0.0,
		 1.0, -1.0, 0.0,
		-1.0,  0.0, 0.0,
		 0.0,  0.0, 0.0,
		 1.0,  0.0, 0.0,
		-1.0,  1.0, 0.0,
		 0.0,  1.0, 0.0,
		 1.0,  1.0, 0.0;

	// element nodal indices
	uMatrix elements(5, 1+4);
	elements <<
		NiHu::quad_1_elem::id, 0, 1, 4, 3,
		NiHu::quad_1_elem::id, 1, 2, 5, 4,
		NiHu::quad_1_elem::id, 3, 4, 7, 6,
		NiHu::tria_1_elem::id, 4, 5, 8, 0,
		NiHu::tria_1_elem::id, 4, 8, 7, 0;

	// create the mesh
	auto msh = NiHu::create_mesh(nodes, elements, NiHu::tria_1_tag(), NiHu::quad_1_tag());

	// create a piecewise constant function space
	auto const &w = NiHu::constant_view(msh);

	// compute number of DOF and allocate result matrix
	size_t nDOF = w.get_num_dofs();
	dMatrix I(nDOF, nDOF);
	I.setZero();

	// create integral operator from kernel and perform weighted double integral
	auto K = NiHu::create_integral_operator(NiHu::laplace_3d_SLP_kernel());
	I <<  w * K[w];

	// Display matrix elements and their sum
	std::cout << "WR matrix:\n" << I << std::endl;
	std::cout << "sum of elements: " << I.sum() << std::endl;

	// Compare to analytical solution
	double anal = 32.0 * (std::log(1.0 + root_two)-(root_two - 1.0)/3.0) / 4.0 / pi;
	std::cout << "log10 error = " << std::log10(std::abs(I.sum() / anal - 1.0)) << std::endl;

	return 0;
}

