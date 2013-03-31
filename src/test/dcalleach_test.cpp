#include "../tmp/control.hpp"

#include <iostream>

class A
{
public:
	A() {  std::cout << "constr" << std::endl; }
	A (A const &other) { std::cout << "copy" << std::endl; }
};

template <class C1, class C2>
struct functor { struct type {
	void operator() (A const &a, A const &b)
	{
	}
};};

int main(void)
{
	A a;
	tmp::d_call_each<
		tmp::vector<A, A, A>,
		tmp::vector<A, A, A>,
		functor<tmp::_1, tmp::_2>,
		A &,
		A &
	>(a, a);
	
	return 0;
}

