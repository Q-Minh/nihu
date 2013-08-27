#include "tmp/integer.hpp"
#include "tmp/vector.hpp"
#include "tmp/print.hpp"

#include <iostream>

int main(void)
{
	typedef tmp::vector<> v0;
	typedef tmp::vector<tmp::int_<1>, tmp::int_<-3>, tmp::int_<0> > v1;
	typedef tmp::vector<tmp::int_<-4>, tmp::int_<26>, tmp::int_<7>, tmp::int_<13> > v2;
	
	std::cout << "Test vectors" << std::endl;
	std::cout << "============" << std::endl;
	std::cout << "v0 : "; print<v0>::eval();
	std::cout << "v1 : "; print<v1>::eval();
	std::cout << "v2 : "; print<v2>::eval();
	std::cout << std::endl;
	
	std::cout << "Testing at" << std::endl;
	std::cout << "==========" << std::endl;
	std::cout << "at<v1, 1> = " << tmp::at<v1, tmp::int_<1> >::type::value << std::endl;
	std::cout << "at<v2, 3> = " << tmp::at<v2, tmp::int_<3> >::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing pop and push" << std::endl;
	std::cout << "====================" << std::endl;
	std::cout << "push_back(v0, 16)   = "; print<tmp::push_back<v0, tmp::int_<16> >::type>::eval();
	std::cout << "push_front(v2, -12) = "; print<tmp::push_front<v2, tmp::int_<-12> >::type>::eval();
	std::cout << "pop_front(v1)       = "; print<tmp::pop_front<v1>::type>::eval();
	std::cout << std::endl;
	
	std::cout << "Testing size" << std::endl;
	std::cout << "============" << std::endl;
	std::cout << "size(v0) = " << tmp::size<v0>::value << std::endl;
	std::cout << std::endl;
	
	// Test vector iterators
	std::cout << "Testing vector iterators" << std::endl;
	std::cout << "========================" << std::endl;
	std::cout << "*begin(v1)       = " << tmp::deref<tmp::begin<v1>::type>::type::value << std::endl;
	std::cout << "*next(begin(v1)) = " << tmp::deref<tmp::next<tmp::begin<v1>::type>::type>::type::value << std::endl;
	std::cout << "*prev(end(v1))   = " << tmp::deref<tmp::prev<tmp::end<v1>::type>::type>::type::value << std::endl;
	std::cout << std::endl;
	
	return 0;
}