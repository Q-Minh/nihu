#ifndef POOL_PATTERN_HPP
#define POOL_PATTERN_HPP

template <class C, unsigned MaxOrder>
class pool
{
public:
	pool(void)
	{
		for (unsigned i = 0; i <= MaxOrder; ++i)
			m_p_data[i] = new C(i);
	}

	~pool(void)
	{
		for (unsigned i = 0; i <= MaxOrder; ++i)
			delete m_p_data[i];
	}

	C const &operator[](unsigned idx) const
	{
		return *m_p_data[idx];
	}

private:
	C *m_p_data[MaxOrder+1];
};

#endif
