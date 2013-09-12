// This file is a part of NiHu, a C++ BEM template library.
// 
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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

#include "util/collection.hpp"
#include <iostream>

class empty
{
};

class A
{
public:
	A(int data) : m_data(data) { }
	int get_data(void) const { return m_data; }
private:
	int m_data;
};

class B
{
public:
	B(double data) : m_data(data) { }
	double get_data(void) const { return m_data; }
private:
	double m_data;
};

int main(void)
{
	typedef collect<empty> empty_collection_t;
	typedef collect<A, empty> A_collection_t;
	typedef collect<B> B_collection_t;

	empty_collection_t ec;
	A_collection_t ac(A(3), empty());
	B_collection_t bc(B(2.0));

	std::cout << ac.A::get_data() << std::endl;

	auto m = merge_data(ac, bc);
	std::cout << m.A::get_data() << std::endl;
}
