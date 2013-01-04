#include "placeholder.hpp"

#include "lambda.hpp"
#include "integer.hpp"

#include <iostream>

int main(void)
{
	std::cout << isPlaceholder<_1>::value << std::endl; // true
	std::cout << isPlaceholder<char>::value << std::endl; // false

	std::cout << std::endl; // false

	std::cout << isPlaceholderExpression< int_<1> >::type::value << std::endl; // false
	std::cout << isPlaceholderExpression< _1 >::value << std::endl; // true
	std::cout << isPlaceholderExpression< next<int_<1> > >::value << std::endl; // false
	std::cout << isPlaceholderExpression< next<_2> >::value << std::endl; // true
	std::cout << isPlaceholderExpression< next< next<_1> > >::value << std::endl; // true
	std::cout << isPlaceholderExpression< plus<int_<2>, int_<3> > >::value << std::endl; // false
	std::cout << isPlaceholderExpression< plus<_1, int_<5> > >::value << std::endl; // true
	std::cout << isPlaceholderExpression< plus<int_<4>, next<_2> > >::value << std::endl; // true
	std::cout << isPlaceholderExpression< plus<_1, _2> >::value << std::endl; // true

	typedef apply< mul<plus<_1,_2>, minus<_1,_2> >, int_<10>, int_<20> >::type res;

	std::cout << std::endl << res::value << std::endl;

	return 0;
}

