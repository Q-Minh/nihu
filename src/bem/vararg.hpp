#include "../tmp/vector.hpp"

template <class T0, class T1, class T2>
class vararg
{
public:
	template <size_t idx>
	struct at_t : tmp::at<tmp::vector<T0, T1, T2>, tmp::int_<idx> > {};

	vararg() {}
	vararg(T0 const &t0, T1 const &t1, T2 const &t2) : m_t0(t0), m_t1(t1), m_t2(t2) {}

	template <size_t idx>
	void set(typename at_t<idx>::type const &);
	template <> void set<0>(typename at_t<0>::type const &t1) { m_t0 = t1; }
	template <> void set<1>(typename at_t<1>::type const &t2) { m_t1 = t2; }
	template <> void set<2>(typename at_t<2>::type const &t3) { m_t2 = t3; }

	template <size_t idx>
	typename at_t<idx>::type const &get(void) const;
	template<> typename at_t<0>::type const &get<0>(void) const { return m_t0; };
	template<> typename at_t<1>::type const &get<1>(void) const { return m_t1; };
	template<> typename at_t<2>::type const &get<2>(void) const { return m_t2; };

private:
	T0 m_t0;
	T1 m_t1;
	T2 m_t2;
};

template<class vararg_t>

