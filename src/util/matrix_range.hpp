#ifndef MATRIX_RANGE_HPP_INCLUDED
#define MATRIX_RANGE_HPP_INCLUDED

namespace iteration
{
	struct plain {};
	struct diagonal {};
}

template <class IterationMode, class OutIt, class InIt>
class matrix_range_impl
{
public:
	class matrix_iterator
	{
	private:
		void preincrement(iteration::plain)
		{
			++m_init;
			if (m_init == m_range.m_iend)
			{
				m_init = m_range.m_ibegin;
				++m_outit;
			}
		}

		void preincrement(iteration::diagonal)
		{
			++m_init;
			++m_outit;
		}

	public:
		matrix_iterator(OutIt outit, InIt init, matrix_range_impl const &range) :
			m_outit(outit), m_init(init), m_range(range)
		{
		}

		matrix_iterator const &operator++()
		{
			preincrement(IterationMode());
			return *this;
		}

		bool operator!=(matrix_iterator const &other) const
		{
			return m_init != other.m_init || m_outit != other.m_outit;
		}

		InIt const &inner(void) const
		{
			return m_init;
		}

		OutIt const &outer(void) const
		{
			return m_outit;
		}

	private:
		OutIt m_outit;
		InIt m_init;
		matrix_range_impl const &m_range;
	};

	matrix_range_impl(OutIt obegin, OutIt oend, InIt ibegin, InIt iend) :
		m_obegin(obegin), m_oend(oend), m_ibegin(ibegin), m_iend(iend)
	{
	}

	matrix_iterator begin(void) const
	{
		return matrix_iterator(m_obegin, m_ibegin, *this);
	}

	matrix_iterator end(void) const
	{
		return end_impl(IterationMode());
	}

private:
	matrix_iterator end_impl(iteration::plain) const
	{
		return matrix_iterator(m_oend, m_ibegin, *this);
	}

	matrix_iterator end_impl(iteration::diagonal) const
	{
		return matrix_iterator(m_oend, m_iend, *this);
	}

	OutIt m_obegin;
	OutIt m_oend;
	InIt m_ibegin;
	InIt m_iend;
};

template <class IterationMode, class OutIt, class InIt>
matrix_range_impl<IterationMode, InIt, OutIt>
	matrix_range(IterationMode, OutIt obegin, OutIt oend, InIt ibegin, InIt iend)
{
	return matrix_range_impl<IterationMode, OutIt, InIt>(obegin, oend, ibegin, iend);
}

#endif // MATRIX_RANGE_HPP_INCLUDED

