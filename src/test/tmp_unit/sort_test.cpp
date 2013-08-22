#include <iostream>
#include "tmp/algorithm.hpp"
#include "tmp/integer.hpp"

template <class Seq, int idx = 0 >
struct print
{
	static void eval(void)
	{
		std::cout << tmp::at<Seq, tmp::int_<idx> >::type::value << ' ';
		print<Seq, idx+1>::eval();
	}
};

template <class Seq>
struct print<Seq, tmp::size<Seq>::type::value>
{
	static void eval(void)
	{
		std::cout << std::endl;
	}
};

int main(void)
{
	typedef tmp::vector<tmp::int_<1>, tmp::int_<3>, tmp::int_<0>, tmp::int_<-2> > vec;
	print<vec>::eval();
	std::cout << "descending: ";
	print<tmp::bubble_sort<vec, tmp::less<tmp::_1, tmp::_2> >::type>::eval();
	std::cout << "ascending: ";
	print<tmp::bubble_sort<vec, tmp::less<tmp::_2, tmp::_1> >::type>::eval();
}
