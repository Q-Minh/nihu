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

#include "library/lib_element.hpp"
#include "core/match_types.hpp"
#include "core/field.hpp"
#include "tmp/control.hpp"
#include <iostream>

typedef NiHu::tria_1_elem test_elem_t;
typedef NiHu::tria_1_elem trial_elem_t;

typedef NiHu::field_view<test_elem_t, NiHu::field_option::isoparametric, NiHu::_1d> test_field_t;
typedef NiHu::field_view<trial_elem_t, NiHu::field_option::isoparametric, NiHu::_1d> trial_field_t;

typedef NiHu::match_type_vector<test_field_t, trial_field_t>::type match_vector_t;

template <class T> struct printer { struct type {
	void operator()(void) { std::cout << T::value << ' '; }
}; };

int main(void)
{
	std::cout << tmp::size<match_vector_t>::value << ":\t";
	tmp::call_each<match_vector_t, printer<tmp::_1>	>();
	std::cout << '\n';
}

