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

#include <iostream>
using namespace std;

template <class Derived>
class crtp_base
{
public:
	Derived const &derived() const
	{
		return static_cast<Derived const &>(*this);
	}
};

template <class T>
class Aimpl
{
public:
	void implementation(void) const
	{
		cout << "a\n";
	}
};

template <class T>
class A :
	public Aimpl<T>,
	public crtp_base<A<T> >
{
};

class B :
	public crtp_base<B>,
	public Aimpl<char>
{
public:
	void implementation(void) const
	{
		cout << "b\n";
	}
};

template <class T>
void fun(crtp_base<T> const &arg)
{
	arg.derived().implementation();
}


int main(void)
{
	A<char> a;
	B const &bref = static_cast<B const &>(static_cast<Aimpl<char> const &>(a));

	a.implementation();
	bref.implementation();

	fun(a);
	fun(bref);

	return 0;
}

