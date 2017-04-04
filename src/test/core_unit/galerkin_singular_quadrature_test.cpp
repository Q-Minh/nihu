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
#include "core/singular_galerkin_quadrature.hpp"

int main(void)
{
	std::cout << "Testing singular Gaussian quadratures" << std::endl << std::endl;
	
	typedef NiHu::singular_galerkin_quadrature<
		NiHu::gauss_family_tag,
		NiHu::line_domain, NiHu::line_domain
		> generator_t;
	generator_t::quadrature_t trial_quad, test_quad;
	generator_t::template generate<NiHu::match::match_0d_type>(test_quad, trial_quad, 1);
	
	for (unsigned i = 0; i < trial_quad.size(); ++i)
	{
		std::cout
			<< trial_quad[i].get_xi() << '\t'
			<< test_quad[i].get_xi() << '\t'
			<< test_quad[i].get_w() << '\t'
			<< std::endl;
	}
	
	return 0;
}

