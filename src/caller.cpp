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
struct Integrate { struct type { void operator () (void) { C::integrate(); } }; };

int main(void)
{
	typedef tiny<Quad, Tria> classes;	// sequence
	
	for_each<classes, Integrate<_1> >();
	
	return 0;
}

