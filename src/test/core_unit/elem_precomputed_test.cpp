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

#include "tmp/sequence.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"
#include "core/element.hpp"
#include "library/lib_element.hpp"

#include <chrono>
#include <vector>

template <class ElemType>
struct tester
{
	struct type
	{
		void operator()(void)
		{
			std::cout << "\n" << NiHu::element_traits::name<ElemType>::value << " (" << ElemType::id << ")\n"
				<< "==========" << std::endl;
			const unsigned  N = 1e6;
			std::vector<ElemType> elements;
			elements.reserve(N);
			// construct elements

			typename ElemType::xi_t xi = ElemType::domain_t::get_center();

			auto tic = std::chrono::system_clock::now();
			for (unsigned i = 0; i < N; ++i)
				elements.push_back(ElemType(ElemType::coords_t::Random()));
			auto toc = std::chrono::system_clock::now();
			auto const_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic).count();


			tic = std::chrono::system_clock::now();
			double n = 0.0;
			for (auto it = elements.begin(); it != elements.end(); ++it)
				n += it->get_normal(xi).norm();
			toc = std::chrono::system_clock::now();
			auto comp_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic).count();

			std::cout << "construction time: " << const_elapsed << std::endl;
			std::cout << "computation time: " << comp_elapsed << std::endl;
			std::cout << "total time: " << const_elapsed + comp_elapsed << std::endl;
		}
	};
};

int main(void)
{
	tmp::call_each<
		tmp::vector<
			NiHu::line_1_elem, NiHu::line_2_elem,
			NiHu::tria_1_elem, NiHu::quad_1_elem,
			NiHu::tria_2_elem, NiHu::quad_2_elem, NiHu::quad_28_elem
		>,
		tester<tmp::_1>
	>();
	return 0;
}
