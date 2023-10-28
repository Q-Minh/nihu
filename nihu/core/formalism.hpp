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

/**
 * \file formalism.hpp
 * \ingroup assembly
 * \brief return weak form formalism from the test and trial field types
 */
#ifndef FORMALISM_HPP
#define FORMALISM_HPP

#include "field.hpp"
#include <type_traits>

namespace NiHu
{

/** \brief definition of different weak form formalisms */
namespace formalism
{
	/** \brief general case when the test field is not Dirac */
	struct general { typedef general type; };
	/** \brief collocational case when the test field is Dirac */
	struct collocational { typedef collocational type; };
}

/** \brief return formalism from Test and Trial field types
 * \tparam TestField the test field type
 * \tparam TrialField the trial field type
 */
template <class TestField, class TrialField, class = void>
struct get_formalism;

/** \brief specialiastion of NiHu::get_formalism for the collocational case */
template <class TestField, class TrialField>
struct get_formalism<TestField, TrialField,	typename std::enable_if<
		field_traits::is_dirac<TestField>::value && !field_traits::is_dirac<TrialField>::value
	>::type
> : formalism::collocational {};


/** \brief specialiastion of NiHu::get_formalism for the general case */
template <class TestField, class TrialField>
struct get_formalism<TestField, TrialField,
	typename std::enable_if<
		!field_traits::is_dirac<TestField>::value && !field_traits::is_dirac<TrialField>::value
	>::type
> : formalism::general {};



/** \brief metafunction to determine if formalism is collocational */
template <class TestField, class TrialField>
struct is_collocational : std::is_same<
	typename get_formalism<TestField, TrialField>::type,
	formalism::collocational
> {};

} // end of namespace NiHu


#endif
