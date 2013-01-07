#include "../tmp/sequence.hpp"

#include <string>
#include <vector>
#include <algorithm>

class empty {};

template <class Begin, class End, class Aggr = empty>
struct generate
{
	typedef typename deref<Begin>::type A;
	struct temp : public A, Aggr
	{
		
	};
 	struct type : public generate<typename next<Begin>::type, End, temp>::type {};
};

template <class End, class Aggr>
struct generate<End, End, Aggr>
{
	struct type : public Aggr {};
};

class A {
public: void f(void) {}
};

class B {
public: void f(void) {}
};

class C {
public: void f(void) {}
};

typedef tiny<std::vector<A>, std::vector<B> > types;

typedef generate<
	typename begin<types>::type,
	typename end<types>::type
>::type Vector;

int main(void)
{
	Vector v;
	v.std::vector<A>::push_back(A());
	v.std::vector<B>::push_back(B());
	
	std::for_each(
		v.std::vector<A>::begin(),
		v.std::vector<A>::end(),
		[] (A &a) { a.f(); }
	);
	
	std::for_each(
		v.std::vector<B>::begin(),
		v.std::vector<B>::end(),
		[] (B &b) { b.f(); }
	);
	
	return 0;
}

