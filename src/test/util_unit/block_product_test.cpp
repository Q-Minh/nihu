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

#include "util/block_product.hpp"
#include <iostream>

int main(void)
{
	// testing a real block product like
	Eigen::Matrix<double, 2, 1> v1 = Eigen::Matrix<double, 2, 1>::Constant(1.0);
	Eigen::Matrix<double, 3, 3> m = Eigen::Matrix<double, 3, 3>::Constant(2.0);
	Eigen::Matrix<double, 4, 1> v2 = Eigen::Matrix<double, 4, 1>::Constant(3.0);;

	// v1 * m * v2^T
	std::cout << block_product(v1, m, v2) << std::endl;

	// testing a pseudo block product
	Eigen::Matrix<double, 1, 1> vv1 = Eigen::Matrix<double, 1, 1>::Constant(1.0);
	double mat = 2.0;
	Eigen::Matrix<double, 1, 1> vv2 = Eigen::Matrix<double, 1, 1>::Constant(3.0);

	std::cout << block_product(vv1, mat, vv2) << std::endl;
}

