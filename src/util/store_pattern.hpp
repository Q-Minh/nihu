#ifndef STORE_PATTERN_HPP
#define STORE_PATTERN_HPP

template <class C>
struct store
{
	static const C m_data;
};

template <class C>
const C store<C>::m_data;

#endif // STORE_PATTERN_HPP
