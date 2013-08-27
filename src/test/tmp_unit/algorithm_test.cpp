#include "tmp/algorithm.hpp"
#include "tmp/integer.hpp"
#include "tmp/vector.hpp"
#include "tmp/print.hpp"

#include <iostream>

int main(void)
{
	typedef tmp::vector<
		tmp::int_<1>,
		tmp::int_<3>,
		tmp::int_<28>,
		tmp::int_<-32>,
		tmp::int_<3>,
		tmp::int_<0>,
		tmp::int_<-2>
	> v1;
	
	typedef tmp::vector<
		tmp::int_<-8>,
		tmp::int_<6>,
		tmp::int_<2>,
		tmp::int_<-1>
	> v2;
	
	
	std::cout << "Testing accumumlation" << std::endl;
	std::cout << "=====================" << std::endl;
	std::cout << "v1       = "; print<v1>::eval();
	std::cout << "v2       = "; print<v2>::eval();
	std::cout << "sum(v1)  = "; std::cout << tmp::accumulate<v1, tmp::int_<0>>::type::value << std::endl;
	std::cout << "prod(v2) = "; std::cout << tmp::accumulate<v2, tmp::int_<1>, tmp::mul<tmp::_1, tmp::_2> >::type::value << std::endl;
	std::cout << "max(v1)  = "; std::cout << tmp::max<v1>::type::value << std::endl;
	std::cout << "min(v1)  = "; std::cout << tmp::min<v1>::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing bubble_sort" << std::endl;
	std::cout << "===================" << std::endl;
	std::cout << "v1   = "; print<v1>::eval();
	std::cout << "v1 < = "; print<tmp::bubble_sort<v1>::type>::eval();
	std::cout << "v1 > = ";print<tmp::bubble_sort<v1, tmp::less<tmp::_1, tmp::_2> >::type>::eval();
	std::cout << std::endl;
	
	std::cout << "Testing concatenate" << std::endl;
	std::cout << "===================" << std::endl;
	std::cout << "v1       = "; print<v1>::eval();
	std::cout << "v2       = "; print<v2>::eval();
	std::cout << "[v2, v1] = "; print<tmp::concatenate<v2, v1>::type>::eval();
	
	std::cout << "Testing find" << std::endl;
	std::cout << "============" << std::endl; 
	std::cout << "v1            = "; print<v1>::eval();
	std::cout << "find(v1, 3)   = "; std::cout << tmp::deref<tmp::find<v1, tmp::int_<3> >::type>::type::value << std::endl;
	std::cout << "find(v1, -32) = "; std::cout << tmp::deref<tmp::find<v1, tmp::int_<-32> >::type>::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << std::boolalpha;
	std::cout << "Testing is_member" << std::endl;
	std::cout << "=================" << std::endl;
	std::cout << "v1          = "; print<v1>::eval();
	std::cout << "{3} in v1 ? : "; 
	std::cout << tmp::is_member<v1, tmp::int_<3> >::type::value << std::endl;
	std::cout << "{2} in v1 ? : "; 
	std::cout << tmp::is_member<v1, tmp::int_<2> >::type::value << std::endl;
	std::cout << "{0} in v1 ? : "; 
	std::cout << tmp::is_member<v1, tmp::int_<0> >::type::value << std::endl;
	std::cout << std::endl;
	
	std::cout << "Testing serialise" << std::endl;
	std::cout << "=================" << std::endl;
	std::cout << "v1           = "; print<v1>::eval();
	std::cout << "v2           = "; print<v2>::eval();
	std::cout << "[v1, v2, v2] = "; print<tmp::serialise<tmp::vector<v1, v2, v2> >::type>::eval();
	std::cout << std::endl;
	
	std::cout << "Testing unique" << std::endl;
	std::cout << "==============" << std::endl;
	std::cout << "v1         = "; print<v1>::eval();
	std::cout << "unique{v1} = "; print<tmp::unique<v1>::type>::eval();
	std::cout << std::endl;
	
	return 0;
}
