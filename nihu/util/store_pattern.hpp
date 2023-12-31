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
 * \file store_pattern.hpp
 * \brief Store static template specialisations
 * \ingroup util
 */
#ifndef STORE_PATTERN_HPP_INCLUDED
#define STORE_PATTERN_HPP_INCLUDED

namespace NiHu
{

/** 
 * \brief Storage class with a static member
 */
template <class C>
struct store
{
    /** 
	 * \brief Return reference to stored data 
	 */
    static C const &get_data(void)
    {
        static const C m_data;
        return m_data;
    }
};

} // end of namespace NiHu

#endif /* STORE_PATTERN_HPP_INCLUDED */
