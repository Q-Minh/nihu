#ifndef PRINT_HPP_INCLUDED
#define PRINT_HPP_INCLUDED

#include <iostream>

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

#endif // PRINT_HPP_INCLUDED
