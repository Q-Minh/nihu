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
 * \file casted_iterator.hpp
 * \brief implementation of an iterator returning static casted values
 */
#ifndef CASTED_ITERATOR_HPP_INCLUDED
#define CASTED_ITERATOR_HPP_INCLUDED

#include <cstddef> // size_t

namespace NiHu
{

/** \brief iterator class provides access to its value_t after static cast
 * \tparam FromIt the original iterator type
 * \tparam To the new value type
 * \tparam Through an intermediate value type if needed for static cast
 */
template <class FromIt, class To, class Through = To>
class casted_iterator :
	private FromIt
{
public:
	/** \brief self returning metafunction */
	typedef casted_iterator type;

	typedef typename FromIt::difference_type difference_type;

	/** \brief the new value type */
	typedef To value_type;
	
	typedef To value_t;
	
	FromIt const &base() const
	{
		return *static_cast<FromIt const *>(this);
	}

	FromIt &base()
	{
		return static_cast<FromIt &>(*this);
	}

	/** \brief (copy) constructor from base iterator
	* \param [in] base base iterator
	*/
	casted_iterator(FromIt const &base = FromIt()) :
		FromIt(base)
	{
	}
	
	bool operator != (casted_iterator const &other) const
	{
		return base() != other.base();
		
	}
	
	bool operator == (casted_iterator const &other) const
	{
		return base() == other.base();
	}
	
	casted_iterator operator+(difference_type offset) const
	{
		return base() + offset;
	}
	
	casted_iterator operator-(difference_type offset) const
	{
		return base() - offset;
	}
	
	difference_type operator-(casted_iterator const &other) const
	{
		return base() - other.base();
	}
	
	casted_iterator &operator++()
	{
		return static_cast<casted_iterator &>(++base());
	}

	casted_iterator operator++(int)
	{
		return base()++;
	}

	/** \brief dereference operator converts dereferenced element to casted type
	* \return the referred casted type class
	*/
	value_t const &operator *(void) const
	{
		return static_cast<value_t const &>(
			static_cast<Through const &>(
				FromIt::operator*()));
	}

	/** \brief index operator converts indexed element to casted type
	* \return the referred casted type class
	*/
	value_t const &operator[] (size_t idx) const
	{
		return static_cast<value_t const &>(
			static_cast<Through const &>(
				base()[idx]));
	}

	/** \brief dereference operator converts dereferenced element to casted type
	* \return the referred casted type pointer
	*/
	value_t const *operator->(void) const
	{
		return static_cast<value_t const*>(
			static_cast<Through const*>(
				&(*base())));
	}
};

}

#endif // CASTED_ITERATOR_HPP_INCLUDED

