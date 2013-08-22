#include <iostream>
#include "tmp/algorithm.hpp"
#include "tmp/integer.hpp"
#include "tmp/vector.hpp"

template <class Seq>
struct print
{
	static void eval(void)
	{
		std::cout << tmp::deref<typename tmp::begin<Seq>::type>::type::value << ", ";
		print<typename tmp::pop_front<Seq>::type>::eval();
	}
};

template <>
struct print<tmp::vector<> >
{
	static void eval(void)
	{
		std::cout << std::endl;
	}
};

int main(void)
{
	typedef tmp::vector<
		tmp::int_<1>,
		tmp::int_<3>,
		tmp::int_<28>,
		tmp::int_<-32>,
		tmp::int_<0>,
		tmp::int_<-2>
	> vec;

	std::cout << "original:\t";
	print<vec>::eval();

	std::cout << "ascending:\t";
	print<tmp::bubble_sort<vec>::type>::eval();

	std::cout << "descending:\t";
	print<tmp::bubble_sort<vec, tmp::less<tmp::_1, tmp::_2> >::type>::eval();
}
