#include <iostream>
#include <typeinfo>
#include "util/store_pattern.hpp"

//! [cache]
template <class T>
class cache
{
public:
	cache() {
		std::cout << "cache constructor " << typeid(T).name() << "\n";
		m_ptr = new T[1000];
		for (size_t i = 0; i < 1000; ++i) m_ptr[i] = i;
	}

	~cache() {
		std::cout << "cache destructor " << typeid(T).name() << "\n";
		delete [] m_ptr;
	}

	T const &operator[](size_t idx) const {
		return m_ptr[idx];
	}

private:
	T *m_ptr;
};
//! [cache]

//! [main]
int main(void)
{
	std::cout << store<cache<int> >::m_data[5] << std::endl;
	std::cout << store<cache<int> >::m_data[25] << std::endl;
	std::cout << store<cache<char> >::m_data[33] << std::endl;
	return 0;
}
//! [main]
