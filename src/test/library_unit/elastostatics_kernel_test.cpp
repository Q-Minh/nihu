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

int main(void)
{
    elastostatics_3d_U_kernel U(.33);

    quad_1_elem::coords_t coords1, coords2;
    coords1 <<
		0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;
    coords2 <<
		2.0, 3.0, 3.0, 1.0,
		0.0, 0.0, 1.0, 1.0,
		0.0, 0.0, 0.0, 0.0;

    quad_1_elem elem1(coords1), elem2(coords2);

    elastostatics_3d_U_kernel::test_input_t test_input(elem1, quad_1_elem::xi_t::Zero());
    elastostatics_3d_U_kernel::trial_input_t trial_input(elem2, quad_1_elem::xi_t::Zero());

    std::cout << U(test_input, trial_input) << std::endl;
}
