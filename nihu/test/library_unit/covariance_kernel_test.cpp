// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2017  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2017  Peter Rucz <rucz@hit.bme.hu>
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

#include "nihu/library/covariance_kernel.hpp"
#include "nihu/library/lib_element.hpp"
#include "nihu/core/double_integral.hpp"

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

template <class K, class E, class Dimension = typename K::dimension_t >
typename NiHu::double_integral<K, NiHu::dirac_field<NiHu::field_view<E, NiHu::field_option::constant, Dimension> >, NiHu::field_view<E, NiHu::field_option::constant, Dimension> >::result_t
/*Eigen::Matrix<typename K::result_t, 1, 1>*/
singular_integral_test(NiHu::kernel_base<K> const &k, E const &e)
{
	// cast the element into a constant field view
	typedef NiHu::field_view<E, NiHu::field_option::constant> F;
	F const &field = static_cast<F const &>(e);

	// compute collocational singular integral
	typedef NiHu::double_integral<K, NiHu::dirac_field<F>, F> di_type;
	return di_type::eval(k, NiHu::dirac(field), field, std::true_type());
}

int main(void)
{
	// two elements shifted (regular case)
    NiHu::line_1_volume_elem::coords_t coords1, coords2;
    coords1 << 	-1., 1.;
    coords2 << 2., 3.;
    NiHu::line_1_volume_elem elem1(coords1), elem2(coords2);

	// Test unidimensional exponential
	typedef NiHu::exponential_covariance_kernel<NiHu::line_1_volume_elem::space_t, NiHu::field_dimension::_1d> exp_kernel_1d_t;
	exp_kernel_1d_t::field_variance_t variance(2.0);
	double d = 2;
	exp_kernel_1d_t C(variance, d);
	std::cout << "Evaluating C kernel in the center of two different elements..." << std::endl;
    std::cout << tester(C, elem1, elem2) << std::endl;
	std::cout << "Analytic: " << std::endl;
	std::cout << variance * std::exp(-2.5/2.) << std::endl;

	std::cout << "Evaluating collocation singular integral over a constant line element..." << std::endl;
	std::cout << "C singular:\n" << singular_integral_test(C, elem1) << std::endl;
	std::cout << "Analytic: " << std::endl;
	std::cout << 2. * d * variance * (1. - std::exp(-1./2.)) << std::endl;

	// Test multidimensional exponential
	Eigen::Matrix<double, 2, 2> matrix;
	matrix << 2.0, 1.0, 1.0, 4.0;
	NiHu::exponential_covariance_kernel<NiHu::line_1_volume_elem::space_t, NiHu::field_dimension::_2d> C2(matrix, d);
	std::cout << tester(C2, elem1, elem2) << std::endl;
	std::cout << singular_integral_test(C2, elem1) << std::endl;

	// Test multidimensional Gaussian
	NiHu::quad_1_volume_elem::coords_t qcoords;
	qcoords << 0.0, 1.0, 1.0, 0.0,
		0.0, 0.0, 1.0, 1.0;
	NiHu::quad_1_volume_elem qelem(qcoords);
	size_t const fdim = 3;
	size_t const sdim = 2;
	Eigen::Matrix<double, fdim, fdim> FieldVar;
	FieldVar.setIdentity();
	Eigen::Matrix<double, sdim, sdim> SpaceVar;
	SpaceVar.setIdentity();
	NiHu::gaussian_covariance_kernel<NiHu::space<double, sdim>, NiHu::field_dimension::_3d> GK(FieldVar, SpaceVar);

	tester(GK, qelem, qelem);
	
    return 0;
}

