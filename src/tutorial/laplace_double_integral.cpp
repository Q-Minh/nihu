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

//! [Includes]
#include "library/lib_element.hpp"
#include "core/weighted_residual.hpp"
#include "library/laplace_kernel.hpp"
//! [Includes]

//! [Typedefs]
// dynamically resizeable double matrix
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
// dynamically resizeable unsigned matrix
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
//! [Typedefs]

//! [Test]
template <class func_space>
void tester(func_space const &w)
{
	// compute number of DOF and allocate result matrix
	int nDOF = w.get_num_dofs();
	dMatrix I(nDOF, nDOF);
	I.setZero();

	// create integral operator from kernel and perform weighted double integral
	auto K = create_integral_operator(laplace_3d_SLP_kernel());
	I <<  w * K[w];

	// Display matrix elements and their sum
	std::cout << "WR matrix:\n" << I << std::endl;
	std::cout << "sum of elements: " << I.sum() << std::endl;

	// Compare to analytical solution
	double anal = 32.0 * (std::log(1.0+std::sqrt(2.0))-(std::sqrt(2.0)-1.0)/3.0) / 4.0/M_PI;
	std::cout << "log10 error = " << std::log10(std::abs(I.sum() / anal - 1.0)) << std::endl;
}
//! [Test]


int main(void)
{
//! [Mesh]
	// nodal coordinates in 3D, 3 nodes
	dMatrix nodes(4,3);
	nodes <<
		-1.0, -1.0, 0.0,
		 1.0, -1.0, 0.0,
		-1.0,  1.0, 0.0,
		 1.0,  1.0, 0.0;

	// element nodal indices
	uMatrix elements(1, 1+4);
	elements << quad_1_elem::id, 0, 1, 3, 2;

	// create the mesh
	auto msh = create_mesh(nodes, elements, tria_1_tag(), quad_1_tag());
//! [Mesh]

//! [Operators]
	// create a piecewise constant function space and call the tester
	std::cout << "Testing with constant field" << std::endl;
	tester(constant_view(msh));

	// create a piecewise linear function space and call the tester
	std::cout << "Testing with isoparametric field" << std::endl;
	tester(isoparametric_view(msh));
//! [Operators]

	return 0;
}

