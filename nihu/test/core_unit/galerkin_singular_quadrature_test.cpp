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

#include <boost/math/constants/constants.hpp>

#include "nihu/core/gaussian_quadrature.hpp"
#include "nihu/core/singular_galerkin_quadrature.hpp"

struct gunc_straight
{
	double operator()(double x, double y)
	{
		return 1.0 / ((x+1) + (y+1));
	}

	double anal()
	{
		return std::log(16.);
	}
};

struct gunc_rect
{
	double operator()(double x, double y)
	{
		return 1.0 / std::sqrt((x+1)*(x+1) + (y+1)*(y+1));
	}

	double anal()
	{
		return 4. * std::asinh(1.);
	}
};

struct func_straight
{
	double operator()(double x, double y)
	{
		return std::log((x+1) + (y+1));
	}

	double anal()
	{
		return 12. * std::log(2.0) - 6;
	}
};


struct func_rect
{
	double operator()(double x, double y)
	{
		return std::log(std::sqrt((x+1)*(x+1) + (y+1)*(y+1)));
	}

	double anal()
	{
		using namespace boost::math::double_constants;
		return pi + std::log(64) - 6;
	}
};


int main(void)
{
	std::cout << "Testing singular Gaussian quadratures" << std::endl << std::endl;
	
	typedef NiHu::singular_galerkin_quadrature<
		NiHu::gauss_family_tag,
		NiHu::line_domain, NiHu::line_domain
		> generator_t;
	gunc_rect f;
	for (int order = 1; order < 15; ++order)
	{
		generator_t::quadrature_t trial_quad, test_quad;
		generator_t::template generate<NiHu::match::match_0d_type>(test_quad, trial_quad, order);
		
		double I = 0;
		for (unsigned i = 0; i < trial_quad.size(); ++i)
		{
			auto y = trial_quad[i].get_xi()(0);
			auto x = test_quad[i].get_xi()(0);
			auto w = trial_quad[i].get_w() * test_quad[i].get_w();
			
			I += f(x, y) * w;
		}
		
		std::cout << order << '\t' << I << '\t' << f.anal() << '\t' << std::log10(std::abs(I/f.anal() - 1)) << std::endl;
	}
	
	return 0;
}

