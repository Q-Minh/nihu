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
	// testing a real block product
	Eigen::Matrix<double, 3, 3> m;
	m.setConstant(2.0);
	Eigen::Matrix<double, 1, 1> v;

	std::cout << "block product:\n"
		<< NiHu::block_product(Eigen::Matrix<double, 2, 1>::Constant(1.0),
		m,
		Eigen::Matrix<double, 2, 1>::Constant(3.0)) << std::endl;

	std::cout << "semi block product:\n"
		<< NiHu::semi_block_product(m, Eigen::Matrix<double, 2, 1>::Constant(3.0)) << std::endl;

	std::cout << "semi block product:\n"
		<< NiHu::semi_block_product(m, v) << std::endl;

	// testing a pseudo block product
	std::cout << "block product:\n"
		<< NiHu::block_product(Eigen::Matrix<double, 2, 1>::Constant(1.0), 2.0, Eigen::Matrix<double, 2, 1>::Constant(3.0)) << std::endl;

	std::cout << "semi block product:\n"
		<< NiHu::semi_block_product(2.0, Eigen::Matrix<double, 2, 1>::Constant(3.0)) << std::endl;
}

