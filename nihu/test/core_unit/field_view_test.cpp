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

#include <iostream>
#include "nihu/core/field.hpp"
#include "nihu/tmp/sequence.hpp"
#include "nihu/tmp/control.hpp"
#include "nihu/tmp/vector.hpp"
#include "nihu/library/lib_element.hpp"

template <class FieldType>
void eval_tester(NiHu::field_base<FieldType> const &field)
{
	std::cout << "DOFs:\t" << field.get_dofs().transpose() << std::endl;
	std::cout << "N(xi0):\t" << FieldType::nset_t::template eval_shape<0>(FieldType::elem_t::domain_t::get_center()).transpose() << std::endl;
}

template <class ElemType, class option>
struct view_tester
{
	struct type
	{
		void operator() (void)
		{
			std::cout << std::endl << "Elem type id:\t" << ElemType::id << std::endl;

			typename ElemType::nodes_t nodes;
			for (unsigned i = 0; i < ElemType::num_nodes; ++i)
				nodes(i) = i;

			ElemType e(ElemType::coords_t::Random(), 0, nodes);
			std::cout << e.get_coords() << std::endl << e.get_id() << std::endl;

			auto const &fv = create_field_view(e, option());
			std::cout << "DOF:\t" << fv.get_dofs().transpose() << std::endl;

			auto const &dfv = dirac(fv);
			std::cout << "Dirac DOF:\t" << dfv.get_dofs().transpose() << std::endl;
		}
	};
};

int main(void)
{
	typedef tmp::vector<NiHu::tria_1_elem, NiHu::quad_1_elem, NiHu::tria_2_elem, NiHu::quad_2_elem> elem_vector;
	typedef tmp::vector<NiHu::field_option::constant, NiHu::field_option::isoparametric> option_vector;

	tmp::d_call_each<
		elem_vector,
		option_vector,
		view_tester<tmp::_1, tmp::_2>
	>();

	NiHu::quad_1_elem e(NiHu::quad_1_elem::coords_t::Random());

	auto const &ifield = NiHu::isoparametric_view(e);
	eval_tester(ifield);

	auto const &cfield = NiHu::constant_view(e);
	eval_tester(cfield);

	return 0;
}

