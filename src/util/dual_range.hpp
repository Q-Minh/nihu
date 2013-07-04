#ifndef DUAL_RANGE_HPP_INCLUDED
#define DUAL_RANGE_HPP_INCLUDED

namespace iteration
{
	struct plain {};
	struct diagonal {};
}


template <class It1, class It2>
class diagonal_iterator :
	public std::pair<It1, It2>
{
	typedef std::pair<It1, It2> base_t;

public:
	diagonal_iterator(It1 it1, It2 it2) : base_t(it1, it2)
	{
	}

	diagonal_iterator const &operator++()
	{
		++base_t::first;
		++base_t::second;
		return *this;
	}

	bool operator!=(diagonal_iterator const &other) const
	{
		return base_t::first != other.get_first() || base_t::second != other.get_second();
	}

	It1 const &get_first(void) const
	{
		return base_t::first;
	}

	It2 const &get_second(void) const
	{
		return base_t::second;
	}
};


template <class It1, class It2>
class matrix_iterator :
	public diagonal_iterator<It1, It2>
{
	typedef diagonal_iterator<It1, It2> base_t;

public:
	matrix_iterator(It1 it1, It2 it2, It2 m_begin, It2 m_end) :
		base_t(it1, it2), m_begin(m_begin), m_end(m_end)
	{
	}

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
	It2 m_begin;
	It2 m_end;
};


template <class Mode, class OutIt, class InIt>
struct iterator_select;

template <class OutIt, class InIt>
struct iterator_select<iteration::plain, OutIt, InIt>
{
	typedef matrix_iterator<OutIt, InIt> type;
};


template <class OutIt, class InIt>
struct iterator_select<iteration::diagonal, OutIt, InIt>
{
	typedef diagonal_iterator<OutIt, InIt> type;
};


template <class IterationMode, class OutIt, class InIt>
class dual_range
{
public:
	typedef typename iterator_select<IterationMode, OutIt, InIt>::type iterator;

	dual_range(OutIt obegin, OutIt oend, InIt ibegin, InIt iend) :
		m_obegin(obegin), m_oend(oend), m_ibegin(ibegin), m_iend(iend)
	{
	}

	iterator begin(void) const
	{
		return begin_impl(IterationMode());
	}

	iterator end(void) const
	{
		return end_impl(IterationMode());
	}

private:
	iterator begin_impl(iteration::plain) const
	{
		if (m_obegin == m_oend || m_ibegin == m_iend)
			return end_impl(iteration::plain());
		return iterator(m_obegin, m_ibegin, m_ibegin, m_iend);
	}

	iterator begin_impl(iteration::diagonal) const
	{
		return iterator(m_obegin, m_ibegin);
	}

	iterator end_impl(iteration::plain) const
	{
		return iterator(m_oend, m_ibegin, m_ibegin, m_iend);
	}

	iterator end_impl(iteration::diagonal) const
	{
		return iterator(m_oend, m_iend);
	}

	OutIt m_obegin;
	OutIt m_oend;
	InIt m_ibegin;
	InIt m_iend;
};



template <class IterationMode, class OutIt, class InIt>
dual_range<IterationMode, OutIt, InIt>
create_dual_range(IterationMode, OutIt obegin, OutIt oend, InIt ibegin, InIt iend)
{
	return dual_range<IterationMode, OutIt, InIt>(obegin, oend, ibegin, iend);
}

#endif // DUAL_RANGE_HPP_INCLUDED

