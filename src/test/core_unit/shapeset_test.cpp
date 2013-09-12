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

#include "tmp/control.hpp"
#include "core/shapeset.hpp"
#include "tmp/vector.hpp"

template <class shape_set>
struct tester
{
	struct type {
		void operator()(void)
		{
			std::cout << std::endl << "shape set id: " << shape_set::id << std::endl;

			std::cout << "corners: " << std::endl;
			for (auto it = shape_set::corner_begin(); it != shape_set::corner_end(); ++it)
				std::cout << it->transpose() << std::endl;
			std::cout << std::endl;

			std::cout << "shape values in corners: " << std::endl;
			for (auto it = shape_set::corner_begin(); it != shape_set::corner_end(); ++it)
				std::cout << shape_set::eval_shape(*it).transpose() << std::endl;
			std::cout << std::endl;

			std::cout << "shape derivative values in corners" << std::endl;
			for (auto it = shape_set::corner_begin(); it != shape_set::corner_end(); ++it)
				std::cout << shape_set::eval_dshape(*it).transpose() << "  sum: " << shape_set::eval_dshape(*it).sum() << std::endl;
			std::cout << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<
			line_0_shape_set, tria_0_shape_set, quad_0_shape_set, brick_0_shape_set,
			line_1_shape_set, tria_1_shape_set, quad_1_shape_set, brick_1_shape_set,
			parallelogram_shape_set, tria_2_shape_set, quad_2_shape_set
		>,
		tester<tmp::_1>
	>();

	return 0;
}

