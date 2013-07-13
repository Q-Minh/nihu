/** \file dual_range.hpp
* \brief implementation of a dual iterator
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
*/
#ifndef DUAL_RANGE_HPP_INCLUDED
#define DUAL_RANGE_HPP_INCLUDED

#include <utility>

/** \brief dual iteration options */
namespace iteration
{
	/** \brief inner and outer iterators (Descartes) */
	struct plain {};
	/** \brief parallel */
	struct diagonal {};
}


/** \brief two iterators traversing in parallel mode
* \tparam It1 the first iterator type
* \tparam It2 the second iterator type
*/
template <class It1, class It2>
class diagonal_iterator :
	public std::pair<It1, It2>
{
	typedef std::pair<It1, It2> base_t;
	
public:
	/** \brief self-returning metafunction */
	typedef diagonal_iterator type;

	/** \brief constructor
	* \param [in] it1 the first iterator
	* \param [in] it2 the second iterator
	*/
	diagonal_iterator(It1 it1, It2 it2) :
		base_t(it1, it2)
	{
	}

	/** \brief increment operator
	* \return reference to the incremented dual iterator
	*/
	diagonal_iterator const &operator++()
	{
		++base_t::first;
		++base_t::second;
		return *this;
	}

	/** \brief compare two iterators
	* \return true if the iterators are different
	*/
	bool operator!=(diagonal_iterator const &other) const
	{
		return base_t::first != other.get_first() ||
			base_t::second != other.get_second();
	}

	/** \brief return first iterator
	* \return reference to the first iterator
	*/
	It1 const &get_first(void) const
	{
		return base_t::first;
	}

	/** \brief return second iterator
	* \return reference to the second iterator
	*/
	It2 const &get_second(void) const
	{
		return base_t::second;
	}
};


/** \brief two iterators traversing in Descartes mode
* \tparam It1 the first iterator type
* \tparam It2 the second iterator type
*/
template <class It1, class It2>
class matrix_iterator :
	public diagonal_iterator<It1, It2>
{
public:
	/** \brief the base type */
	typedef diagonal_iterator<It1, It2> base_t;

	/** \brief self-returning metafunction */
	typedef matrix_iterator type;

	/** \brief constructor
	* \param [in] it1 the first iterator
	* \param [in] it2 the second iterator
	* \param [in] begin2 begin of the second iterator
	* \param [in] end2 end of the second iterator
	*/
	matrix_iterator(It1 it1, It2 it2, It2 begin2, It2 end2) :
		base_t(it1, it2), m_begin(begin2), m_end(end2)
	{
	}

	/** \brief increment operator
	* \return reference to the incremented dual iterator
	*/
	matrix_iterator const &operator++()
	{
		++base_t::second;
		if (base_t::second == m_end)
		{
			base_t::second = m_begin;
			++base_t::first;
		}
		return *this;
	}

private:
	/** \brief begin of the second iterator */
	It2 m_begin;
	/** \brief end of the second iterator */
	It2 m_end;
};


/** \brief metafunction assigning a dual iterator to an iteration mode */
template <class Mode, class OutIt, class InIt>
struct iterator_select;

/** \brief specialisation of ::iterator_select for the plain case */
template <class OutIt, class InIt>
struct iterator_select<iteration::plain, OutIt, InIt> :
	public matrix_iterator<OutIt, InIt> {};

/** \brief specialisation of ::iterator_select for the diagonal case */
template <class OutIt, class InIt>
struct iterator_select<iteration::diagonal, OutIt, InIt> :
	public diagonal_iterator<OutIt, InIt> {};

/** \brief a combination of two ranges */
template <class IterationMode, class OutIt, class InIt>
class dual_range
{
public:
	/** \brief the range's iteartor type */
	typedef typename iterator_select<IterationMode, OutIt, InIt>::type iterator;

	/** \brief constructor
	* \param [in] obegin outer begin
	* \param [in] oend outer end
	* \param [in] ibegin inner begin
	* \param [in] iend inner end
	*/
	dual_range(OutIt obegin, OutIt oend, InIt ibegin, InIt iend) :
		m_obegin(obegin), m_oend(oend), m_ibegin(ibegin), m_iend(iend)
	{
	}

	/** \brief return begin iterator of dual range */
	iterator begin(void) const
	{
		return begin_impl(IterationMode());
	}

	/** \brief return end iterator of dual range */
	iterator end(void) const
	{
		return end_impl(IterationMode());
	}

private:
	/** \brief begin iterator of the plain range */
	iterator begin_impl(iteration::plain) const
	{
		if (m_obegin == m_oend || m_ibegin == m_iend)
			return end_impl(iteration::plain());
		return iterator(m_obegin, m_ibegin, m_ibegin, m_iend);
	}

	/** \brief begin iterator of the diagonal range */
	iterator begin_impl(iteration::diagonal) const
	{
		return iterator(m_obegin, m_ibegin);
	}

	/** \brief end iterator of the plain range */
	iterator end_impl(iteration::plain) const
	{
		return iterator(m_oend, m_ibegin, m_ibegin, m_iend);
	}

	/** \brief end iterator of the diagonal range */
	iterator end_impl(iteration::diagonal) const
	{
		return iterator(m_oend, m_iend);
	}

	/** \brief the outer begin */
	OutIt m_obegin;
	/** \brief the outer end */
	OutIt m_oend;
	/** \brief the inner begin */
	InIt m_ibegin;
	/** \brief the inner end */
	InIt m_iend;
};



/** \brief factory function to create a dual range
* \tparam IterationMode the iteration mode
* \tparam OutIt the outer iterator type
* \tparam InIt the inner iterator type
* \param [in] obegin the outer iterator's begin
* \param [in] oend the outer iterator's end
* \param [in] ibegin the inner iterator's begin
* \param [in] iend the inner iterator's end
* \return dual range constructed from the iterators
*/
template <class IterationMode, class OutIt, class InIt>
dual_range<IterationMode, OutIt, InIt>
	create_dual_range(IterationMode, OutIt obegin, OutIt oend, InIt ibegin, InIt iend)
{
	return dual_range<IterationMode, OutIt, InIt>(obegin, oend, ibegin, iend);
}

#endif // DUAL_RANGE_HPP_INCLUDED

