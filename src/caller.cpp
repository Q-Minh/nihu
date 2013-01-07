#include "tmp/sequence.hpp"
#include "tmp/call_each.hpp"
#include <iostream>

struct Quad {
	static void integrate (void) { std::cout << "Quad" << std::endl; }
	static void integrate (char c) { std::cout << "Quad " << c << std::endl; }
};
struct Tria {
	static void integrate (void) { std::cout << "Tria" << std::endl; }
	static void integrate (char c) { std::cout << "Tria " << c << std::endl; }
};


template <class C>
struct Integrate { void operator () (void) { C::integrate(); } };

template <class C>
struct Integrate2 { void operator () (char c) { C::integrate(c); } };

int main(void)
{
	typedef tiny<Quad, Tria> classes;	// sequence
	typedef begin<classes>::type begin;	// begin iterator
	typedef end<classes>::type end;		// end iterator
	
	typedef call_each<begin, end, Integrate<_1> >::type integrate_over_classes;
	typedef call_each<begin, end, Integrate2<_1>, char>::type integrate2_over_classes;

	integrate_over_classes::apply();
	integrate2_over_classes::apply('a');

	return 0;
}

