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
 * \file pool_pattern.hpp
 * \brief Store an array of instances
 * \ingroup util
 */
#ifndef POOL_PATTERN_HPP_INCLUDED
#define POOL_PATTERN_HPP_INCLUDED

#include <stdexcept>

namespace NiHu
{

/**
 * \brief Class storing a vector of class instances
 * \tparam C the stored class
 * \tparam MaxOrder the array size
 * \details 
 * In order to create a pool, the stored type \c C must have a constructor that 
 * can take one integer as a parameter.
 */
template <class C, unsigned MaxOrder>
class pool
{
public:
	/** \brief constructor */
	pool(void)
	{
		for (unsigned i = 0; i <= MaxOrder; ++i)
			m_p_data[i] = new C(i);
	}

	/** \brief destructor */
	~pool(void)
	{
		for (unsigned i = 0; i <= MaxOrder; ++i)
			delete m_p_data[i];
	}

	/** 
	 * \brief Index operator
	 * \param [in] idx Element index
	 * \returns The indexed element
	 */
	C const &operator[](unsigned idx) const
	{
		if (idx > MaxOrder)
			throw std::out_of_range("Pool overindexing");
		return *m_p_data[idx];
	}

private:
	C *m_p_data[MaxOrder+1];
};

} // end of namespace NiHu

#endif /* POOL_PATTERN_HPP_INCLUDED */
