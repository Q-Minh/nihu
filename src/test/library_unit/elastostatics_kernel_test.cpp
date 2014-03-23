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

#include "library/elastostatics_kernel.hpp"
#include "library/lib_element.hpp"

template <class K, class E1, class E2>
typename K::result_t tester(K const &k, E1 const &e1, E2 const &e2)
{
	// kernel inputs for U kernel
    typename K::test_input_t test_input(e1, E1::domain_t::get_center());
    typename K::trial_input_t trial_input(e2, E2::domain_t::get_center());
    return k(test_input, trial_input);
}

int main(void)
{
	// two elements shifted (regular case)
    quad_1_elem::coords_t coords1, coords2;
    coords1 <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
    coords2 = coords1;
    coords2.row(1) += Eigen::Matrix<double, 1, 4>::Constant(2.0);

    quad_1_elem elem1(coords1), elem2(coords2);
    
    std::cout << tester(elastostatics_3d_U_kernel(.33), elem1, elem2) << std::endl;
	std::cout << tester(elastostatics_3d_T_kernel(.33), elem1, elem2) << std::endl;
    std::cout << tester(create_couple_kernel(elastostatics_3d_U_kernel(.33), elastostatics_3d_T_kernel(.33)), elem1, elem2) << std::endl;
    
    return 0;
}

