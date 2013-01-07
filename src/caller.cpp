#include "sequence.hpp"
#include "call_each.hpp"
#include <iostream>

struct Quad {
	static void noparam (void) { std::cout << "Quad" << std::endl; }
	static void param (char c) { std::cout << "Quad " << c << std::endl; }
};
struct Tria {
	static void noparam (void) { std::cout << "Tria" << std::endl; }
	static void param (char c) { std::cout << "Tria " << c << std::endl; }
};


template <class C>
struct Noparam { void operator () (void) { C::noparam(); } };

template <class C>
struct Param { void operator () (char c) { C::param(c); } };

int main(void)
{
	typedef tiny<Quad, Tria> classes;
	typedef begin<classes>::type begin;
	typedef end<classes>::type end;
	
	typedef call_each<begin, end, Noparam<_1> >::type noparam_over_classes;
	typedef call_each<begin, end, Param<_1>, char>::type param_over_classes;
	
	noparam_over_classes::apply();
	param_over_classes::apply('a');

	return 0;
}

