#include "../tmp/algorithm.hpp"
#include "../tmp/sequence.hpp"
#include "../tmp/call_each.hpp"

#include <string>
#include <vector>
#include <algorithm>

#include <iostream>

class A { public: void f(void) { std::cout << 'a' <<std::endl; } };
class B { public: void f(void) { std::cout << 'b' <<std::endl; } };
class C { public: void f(void) { std::cout << 'c' <<std::endl; } };

typedef tiny<A, B, C> types;

// convert type sequence into container sequence
typedef transform<
	begin<types>::type,
	end<types>::type,
	inserter<tiny<>, push_back<_1,_2> >,
	vectorize<_1>
>::type vtypes;

// generate big Container from containers
typedef inherit<
	typename begin<vtypes>::type,
	typename end<vtypes>::type
>::type Vector;

template  <class T>
struct F
{
	void operator() (Vector &v)
	{
		std::for_each(
			v.T::begin(),
			v.T::end(),
			[] (typename T::value_type t) { t.f(); }
		);
	}
};

int main(void)
{
	Vector v;

	v.std::vector<C>::push_back(C());
	v.std::vector<A>::push_back(A());
	v.std::vector<B>::push_back(B());
	v.std::vector<A>::push_back(A());

	typedef call_each<
		begin<vtypes>::type,
		end<vtypes>::type,
		F<_1>,
		Vector &
	>::type caller; // Caller csucsu

	caller::apply(v);
	
	return 0;
}

