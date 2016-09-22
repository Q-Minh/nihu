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

#include "library/elastodynamics_kernel.hpp"
#include "library/lib_element.hpp"
#include "core/double_integral.hpp"

template <class K, class E1, class E2>
typename K::result_t tester(K const &k, E1 const &e1, E2 const &e2)
{
	// instantiate a kernel input x1 in the center of elem e1
	typename K::test_input_t x1(e1, E1::domain_t::get_center());
	// instantiate a kernel input x2 in the center of elem e2
    typename K::trial_input_t x2(e2, E2::domain_t::get_center());
	// ealuate the kernel K(x1, x2)
    return k(x1, x2);
}


template <class K, class E>
typename K::result_t singular_integral_test(kernel_base<K> const &k, E const &e)
{
	// cast the element into a constant field view
	typedef field_view<E, field_option::constant> F;
	F const &field = static_cast<F const &>(e);

	// compute collocational singular integral
	typedef double_integral<K, dirac_field<F>, F> di_type;
	return di_type::eval(k, dirac(field), field, std::true_type());
}

int main(void)
{
	// two elements shifted (regular case)
    quad_1_elem::coords_t coords1, coords2;
    coords1 <<
		-1.,  1,  1., -1.,
		-1., -1., 1.,  1.,
		 0.,  0., 0.,  0.;
    coords2 = coords1;
    coords2.row(1) += Eigen::Matrix<double, 1, 4>::Constant(2.0);

    quad_1_elem elem1(coords1), elem2(coords2);
	double nu = 0.3;
	double rho = 100;
	double mu = 1.e8;
	double om = 2.*M_PI*10.;

	elastodynamics_3d_U_kernel U(nu, rho, mu, om);
	elastodynamics_3d_T_kernel T(nu, rho, mu, om);

	std::cout << "Evaluating U kernel in the center of two different elements..." << std::endl;
    std::cout << tester(U, elem1, elem2) << std::endl;
	std::cout << "Evaluating T kernel in the center of two different elements..." << std::endl;
    std::cout << tester(T, elem1, elem2) << std::endl;

	std::cout << "Evaluating collocation singular integral over a constant quad element..." << std::endl;
	std::cout << "U singular:\n" << singular_integral_test(U, elem1) << std::endl;
	std::cout << "Evaluating collocation singular integral over a constant quad element..." << std::endl;
	std::cout << "T singular:\n" << singular_integral_test(T, elem1) << std::endl;

    return 0;
}

