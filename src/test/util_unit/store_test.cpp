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
	std::cout << store<cache<int> >::get_data()[5] << std::endl;
	std::cout << store<cache<int> >::get_data()[25] << std::endl;
	std::cout << store<cache<char> >::get_data()[33] << std::endl;
	return 0;
}
//! [main]

