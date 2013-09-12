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

#include "tmp/lambda.hpp"
#include "tmp/integer.hpp"
#include "tmp/vector.hpp"
#include "tmp/control.hpp"
#include "tmp/print.hpp"

#include <iostream>
#include <typeinfo>

// Metafunction with 1 parameter
template <class A>
struct MF1 
{
	typedef A type;
};

// Metafunction with 2 parameters
template <class A, class B>
struct MF2 {
	typedef B type;
};

// Transformation with 1 parameter
template <class A>
struct TR1 
{ 
	struct type
	{
		void operator() (void)
		{
			std::cout <<  A::value << "; ";
		}
	};
};

// Transformation with 2 parameters
template <class A, class B>
struct TR2 
{ 
	struct type
	{
		void operator() (void)
		{
			std::cout << A::value << " " << B::value << "; ";
		}
	};
};

int main(void)
{
	// Predeclaration of vector types
	// typedef tmp::vector<tmp::int_<1>, tmp::int_<2>, tmp::int_<3> > v1;
	// typedef tmp::vector<tmp::int_<-11>, tmp::int_<-23>, tmp::int_<-38> > v2;
	
	// Argument selector test
	std::cout << "Testing argument selector" << std::endl;
	std::cout << "=========================" << std::endl;
	std::cout << "<1, 0, 2> [1] : " << tmp::select_argument<1, tmp::int_<1>, tmp::int_<0>, tmp::int_<2> >::type::value << std::endl;
	std::cout << "<1, 0, 2> [3] : " << tmp::select_argument<3, tmp::int_<1>, tmp::int_<0>, tmp::int_<2> >::type::value << std::endl;
	std::cout << std::endl;
	
	// isPlaceholder test
	std::cout << std::boolalpha;
	std::cout << "Testing isPlaceholder" << std::endl;
	std::cout << "=====================" << std::endl;
	std::cout << "<int> : " << tmp::isPlaceholder<int>::value << std::endl;
	std::cout << "<_2>  : " << tmp::isPlaceholder<tmp::_2>::value << std::endl;
	std::cout << std::endl;
	
	// isPlaceholderExpression test
	std::cout << "Testing isPlaceholderExpression" << std::endl;
	std::cout << "===============================" << std::endl;
	std::cout << "<bool>                 : " << tmp::isPlaceholderExpression<bool>::value << std::endl;
	std::cout << "<_1>                   : " << tmp::isPlaceholderExpression<tmp::_1>::value << std::endl;
	std::cout << "<MF1<int> >            : " << tmp::isPlaceholderExpression<MF1<int > >::value << std::endl;
	std::cout << "<MF1<_1> >             : " << tmp::isPlaceholderExpression<MF1<tmp::_1> >::value << std::endl;
	std::cout << "<apply<MF1<_1>, int> > : " << tmp::isPlaceholderExpression<tmp::apply<MF1<tmp::_1>, int >::type >::value << std::endl;
	std::cout << std::endl;
	
	// lambda test
	std::cout << "Testing lambda" << std::endl;
	std::cout << "==============" << std::endl;
	
	auto app_t2 = tmp::apply<TR2<tmp::_1, tmp::_2>, tmp::int_<-6>, tmp::int_<8> >::type();
	std::cout << "apply<TR2, -6, 8>  : "; app_t2(); std::cout << std::endl;
	
	auto lam1 = tmp::lambda<TR1<tmp::_1 > >::type::apply<tmp::int_<1> >::type();
	auto lam2 = tmp::lambda<TR2<tmp::_1, tmp::_2> >::type::apply<tmp::int_<-4>, tmp::int_<3> >::type();
	
	std::cout << "lambda<TR1, _1>::apply<1> : "; lam1(); std::cout << std::endl;
	std::cout << "lambda<TR2, _1, _2>::apply<-4, 3> : "; lam2(); std::cout << std::endl;
	std::cout << std::endl;
	
	return 0;
}
