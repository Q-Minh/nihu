#include "../tmp/lambda.hpp"
#include "../tmp/integer.hpp"
#include "../tmp/vector.hpp"
#include "../tmp/control.hpp"
#include <iostream>
#include <typeinfo>

using namespace tmp;

template <class A>
struct MF1 {
	typedef A type;
};

template <class A, class B>
struct MF2 {
	typedef B type;
};

int main(void)
{
	// test isPlaceholder
	std::cout << "Testing tmp::isPlaceholder" << std::endl;
	std::cout << "==========================" << std::endl;
	
	static_assert(isPlaceholder<_1>::type::value == true, "isPlaceholder error");
	std::cout << "isPlaceholder<_1> = " << std::boolalpha << isPlaceholder<_1>::type::value << std::endl;
	
	static_assert(isPlaceholder<int>::type::value == false, "isPlaceholder error");
	std::cout << "isPlaceholder<_1> = " << std::boolalpha << isPlaceholder<int>::type::value << std::endl;
	
	// test isPlaceholderExpression
	std::cout << std::endl;
	std::cout << "Testing tmp::isPlaceholderExpression" << std::endl;
	std::cout << "====================================" << std::endl;
	
	std::cout << "isPlaceholderExpression<_1> = " << std::boolalpha << isPlaceholderExpression<_1>::type::value << std::endl;
	
	std::cout << "isPlaceholderExpression<int> = " << std::boolalpha << isPlaceholderExpression<int>::type::value << std::endl;
	
	// test containsPlaceholderExpression
	std::cout << std::endl;
	std::cout << "Testing tmp::containsPlaceholderExpression" << std::endl;
	std::cout << "==========================================" << std::endl;
	
	std::cout << "<_1> = " << std::boolalpha << containsPlaceholderExpression<_1>::type::value << std::endl;

	std::cout << "<int> = " << std::boolalpha << containsPlaceholderExpression<int>::type::value << std::endl;

	std::cout << "<_1, int> = " << std::boolalpha << containsPlaceholderExpression<_1, int>::type::value << std::endl;

	std::cout << "<MF1< _1 > > = " << std::boolalpha << containsPlaceholderExpression<MF1<_1> >::type::value << std::endl;

	std::cout << "<MF1< int > > = " << std::boolalpha << containsPlaceholderExpression<MF1< int > >::type::value << std::endl;

	std::cout << "<MF2< int, _1 > > = " << std::boolalpha << containsPlaceholderExpression<MF2< int, _1 > >::type::value << std::endl;

	std::cout << "<MF2< int, char > > = " << std::boolalpha << containsPlaceholderExpression<MF2< int, char > >::type::value << std::endl;

	
	return 0;
}