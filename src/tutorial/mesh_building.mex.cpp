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

#include "core/mesh.hpp"

//! [Typedefs]
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic> uMatrix;
//! [Typedefs]

void test1()
{
//! [Matrices]
	dMatrix nodes(6,3);
	nodes <<				//	3   4   5
		0.0, 0.0, 0.0,		//	+---+---+
		1.0, 0.0, 0.0,		//	|   |2 /|
		2.0, 0.0, 0.0,		//	| 0 | / |
		0.0, 1.0, 0.0,		//	|   |/ 1|
		1.0, 1.0, 0.0,		//	+---+---+
		2.0, 1.0, 0.0;		//	0   1   2

	uMatrix elements(3,5);
	elements <<
		quad_1_elem::id, 0, 1, 4, 3,
		tria_1_elem::id, 1, 2, 5, 0,
		tria_1_elem::id, 1, 5, 4, 0;
//! [Matrices]

//! [Creation]
	auto my_mesh = create_mesh(nodes, elements, _quad_1_tag(), _tria_1_tag());
//! [Creation]
}

void test2(void)
{
//! [2D Circle]
	double R = 1.0;				// radius
	unsigned N = 180;			// number of elements

	dMatrix nodes(N, 2);		// nodal locations
	for (unsigned i = 0; i < N; ++i)
	{
		double phi = i * (2*M_PI/N);
		nodes(i,0) = R * cos(phi);
		nodes(i,1) = R * sin(phi);
	}

	uMatrix elements(N,3);		// element connections
	for (unsigned i = 0; i < N; ++i)
	{
		elements(i,0) = line_1_elem::id;
		elements(i,1) = i;
		elements(i,2) = (i+1)%N;
	}

	auto mesh = create_mesh(nodes, elements, _line_1_tag());	// homogeneous mesh
//! [2D Circle]
}



//! [Matlab example]
#include "util/mex_matrix.hpp"

// the MATLAB entry point
void mexFunction(int nlhs, mxArray *lhs[], int nrhs, mxArray const *rhs[])
{
	// matrix import
	mex::real_matrix<double> nodes(rhs[0]), elements(rhs[1]);
	// create the mesh from Matlab
	auto surf_mesh = create_mesh(nodes, elements, _tria_1_tag(), _quad_1_tag());
}
//! [Matlab example]


int main(void)
{
	test1();
	test2();
}

