/**
 * @file elem_center_iterator.hpp
 * @brief Classes for iterating through elements and fields
 * @ingroup fmm_clus
 */

#ifndef ELEM_CENTER_ITERATOR_HPP_INCLUDED
#define ELEM_CENTER_ITERATOR_HPP_INCLUDED

#include <iterator>

namespace NiHu
{
namespace fmm
{

/// \brief iterate through element centers
/// \tparam It the iterator iterating through elements
template <class It>
class elem_center_iterator
	: public std::iterator<std::forward_iterator_tag, typename It::value_type::x_t>
	, private It
{
	typedef typename It::value_type elem_t;
	typedef typename elem_t::x_t x_t;

public:
	/// \brief conversion constructor
	/// \param [in] it the element iterator
	elem_center_iterator(It it = It())
		: It(it)
	{
	}

	/// \brief  helper function to reach the element iterator interface
	/// \return the object as an element iterator
	It const &base() const
	{
		return static_cast<It const &>(*this);
	}

	/// \brief  helper function to reach the element iterator interface
	/// \return the object as an element iterator
	It &base()
	{
		return static_cast<It &>(*this);
	}

	/// \brief difference of two iterators
	/// \param [in] other the other iterator
	/// \return the iterator difference
	ptrdiff_t operator-(elem_center_iterator const &other) const
	{
		return base() - other;
	}

	/// \brief iterator arithmetics
	/// \param [in] diff the offset
	/// \return the offseted iterator
	elem_center_iterator operator-(ptrdiff_t diff) const
	{
		return base() - diff;
	}

	/// \brief iterator arithmetics
	/// \param [in] diff the offset
	/// \return the offseted iterator
	elem_center_iterator operator+(ptrdiff_t diff) const
	{
		return base() + diff;
	}

	/// \brief preincrenet operator
	/// \return reference to the incremented instance
	elem_center_iterator &operator++()
	{
		return static_cast<elem_center_iterator &>(++base());
	}

	/// \brief postincrenet operator
	/// \return the original instance
	elem_center_iterator operator++(int)
	{
		return base()++;
	}

	/// \brief dereference operator
	/// \return the element center
	x_t const &operator*() const
	{
		return (*base()).get_center();
	}

	/// \brief index operator
	/// \return the element center of the offseted iterator
	x_t const &operator[](size_t idx) const
	{
		return base()[idx].get_center();
	}

	/// \brief non-equality operator
	/// \return true if two iterators are different
	bool operator != (elem_center_iterator const& other) const
	{
		return base() != other;
	}

	/// \brief equality operator
	/// \return true if two iterators are different
	bool operator == (elem_center_iterator const& other) const
	{
		return base() == other;
	}
};

/// \brief factory function to create an element center iterator
/// \tparam It the iterator class
/// \param [in] it the iterator instance
/// \return the elem center iterator instance
template <class It>
auto create_elem_center_iterator(It it)
{
	return elem_center_iterator<It>(it);
}


/// \brief iterate through field centers
/// \tparam It the iterator iterating through fields
template <class It>
class field_center_iterator
	: public std::iterator<std::forward_iterator_tag, typename It::value_type::elem_t::x_t>
	, private It
{
	typedef typename It::value_type field_t;
	typedef typename field_t::elem_t elem_t;
	typedef typename elem_t::x_t x_t;

public:
	/// \brief conversion constructor
	/// \param [in] it the element iterator
	field_center_iterator(It it = It())
		: It(it)
	{
	}

	/// \brief  helper function to reach the field iterator interface
	/// \return the object as an field iterator
	It const &base() const
	{
		return *static_cast<It const *>(this);
	}

	/// \brief  helper function to reach the field iterator interface
	/// \return the object as an field iterator
	It &base()
	{
		return *static_cast<It *>(this);
	}

	/// \brief difference of two iterators
	/// \param [in] other the other iterator
	/// \return the iterator difference
	ptrdiff_t operator-(field_center_iterator const &other) const
	{
		return base() - other;
	}

	/// \brief iterator arithmetics
	/// \param [in] diff the offset
	/// \return the offseted iterator
	field_center_iterator operator-(ptrdiff_t diff) const
	{
		return base() - diff;
	}

	/// \brief iterator arithmetics
	/// \param [in] diff the offset
	/// \return the offseted iterator
	field_center_iterator operator+(ptrdiff_t diff) const
	{
		return base() + diff;
	}

	/// \brief preincrenet operator
	/// \return the incremented instance
	field_center_iterator &operator++()
	{
		return static_cast<field_center_iterator &>(++base());
	}

	/// \brief postincrenet operator
	/// \return the original instance
	field_center_iterator operator++(int)
	{
		return base()++;
	}

	/// \brief dereference operator
	/// \return the element center
	x_t const &operator*() const
	{
		return (*base()).get_elem().get_center();
	}

	/// \brief index operator
	/// \return the element center of the offseted iterator
	x_t const &operator[](size_t idx) const
	{
		return base()[idx].get_elem().get_center();
	}

	/// \brief non-equality operator
	/// \return true if two iterators are different
	bool operator != (field_center_iterator const& other) const
	{
		return base() != other;
	}

	/// \brief equality operator
	/// \return true if two iterators are different
	bool operator == (field_center_iterator const& other) const
	{
		return base() == other;
	}
};

/// \brief factory function to create a field center iterator
/// \tparam It the iterator class
/// \param [in] it the iterator instance
/// \return a field center iterator
template <class It>
auto create_field_center_iterator(It it)
{
	return field_center_iterator<It>(it);
}


} // end of namespace fmm
} // end of namespace NiHu


#endif /* ELEM_CENTER_ITERATOR_HPP_INCLUDED */
