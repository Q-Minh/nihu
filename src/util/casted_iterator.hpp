/**
 * \file casted_iterator.hpp
 */
#ifndef CASTED_ITERATOR_HPP_INCLUDED
#define CASTED_ITERATOR_HPP_INCLUDED

/** \brief iterator class provides access to the mesh's elements as fields */
template <class FromIt, class To, class Through = To>
class casted_iterator : public FromIt
{
public:
	/** \brief the pointed data type */
	typedef To value_t;

	/** \brief (copy) constructor from base iterator
	* \param [in] base base iterator
	*/
	casted_iterator(FromIt const &base) :
		FromIt(base)
	{
	}

	/** \brief dereference operator converts dereferenced element into field_view
	* \return the referred field_view class
	*/
	value_t const &operator *(void) const
	{
		return static_cast<value_t const &>(
			static_cast<Through const &>(
				FromIt::operator*()));
	}

	/** \brief pointer dereference operator
	* \return pointer to the referred field view instance
	*/
	value_t const *operator->(void) const
	{
		return &(*(*this));
	}
};

#endif // CASTED_ITERATOR_HPP_INCLUDED

