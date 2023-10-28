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

#include "core/gaussian_quadrature.hpp"
#include "core/duffy_quadrature.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"

template <class Family, class Domain>
struct tester
{
	struct type
	{
		void operator() (void)
		{
			// instantiate a quadrature
			typename NiHu::quadrature_type<Family, Domain>::type q(5);
			// print points and weights
			std::cout << q << std::endl;
		}
	};
};


void test_transform(void)
{
	NiHu::gaussian_quadrature<NiHu::quad_domain> q(4);
	Eigen::Matrix<double, 4, 2> coords;
	coords <<
		0.0, 0.0,
		1.0, 0.0,
		1.0, 1.0,
		1.0, 1.0;
	std::cout << q.transform<NiHu::quad_1_shape_set>(coords) << std::endl;
}


template <class Family, class LSet>
struct duffy_test
{
	struct type {
		void operator()(void)
		{
			typedef NiHu::duffy_quadrature<Family, LSet> duffy_t;

			std::cout << duffy_t::on_corner(3, 2) << std::endl;

			std::cout << duffy_t::on_face(3, LSet::domain_t::get_center()) << std::endl;
		}
	};
};


int main(void)
{
	std::cout << "Testing regular quadratures" << std::endl << std::endl;
	tmp::call_each<
		tmp::vector<NiHu::line_domain, NiHu::tria_domain, NiHu::quad_domain>,
		tester<NiHu::gauss_family_tag, tmp::_1>
	>();

	std::cout << "Testing quadrature transforms" << std::endl << std::endl;
	test_transform();

	std::cout << "Testing Duffy quadratures" << std::endl << std::endl;
	tmp::call_each<
		tmp::vector<NiHu::tria_1_shape_set, NiHu::quad_1_shape_set>,
		duffy_test<NiHu::gauss_family_tag, tmp::_1>
	>();

	return 0;
}

