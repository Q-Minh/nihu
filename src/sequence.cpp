#include "sequence.hpp"

#include "algorithm.hpp"

#include <iostream>

int main(void)
{
	typedef tiny<int_<1>, int_<10>, int_<20> > seq0;
	
	// at
	std::cout <<
		at<seq0, int_<0> >::type::value << " " <<
		at<seq0, int_<1> >::type::value << " " <<
		at<seq0, int_<2> >::type::value << " " <<
		std::endl;
	
	// size
	std::cout << size<seq0>::type::value << std::endl;

	// begin, next, deref
	std::cout <<
		deref<
			next<next<begin<seq0>::type>::type>::type
		>::type
		::value << std::endl;
		
	// end, prev, deref
	std::cout <<
		deref<
			prev<prev<end<seq0>::type>::type>::type
		>::type
		::value << std::endl;
		
	// clear, push_front, push_back
	typedef push_back<
		push_front<
			clear<seq0>::type,
			int_<-1>
		>::type, int_<-2>
	>::type seq1;
	
	std::cout <<
		at<seq1, int_<0> >::type::value << " " <<
		at<seq1, int_<1> >::type::value << " " <<
		std::endl;
		
	// sum by accumulation
	std::cout <<
		accumulate<
			begin<seq0>::type,
			end<seq0>::type,
			int_<0> // default summation
		>::type
		::value << std::endl;
	
	// factorial by accumulation
	std::cout <<
		accumulate<
			begin<seq0>::type,
			end<seq0>::type,
			int_<1>,
			mul<_1,_2>
		>::type
		::value << std::endl;
	
	return 0;
}


