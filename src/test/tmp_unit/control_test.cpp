#include "tmp/algorithm.hpp"
#include "tmp/control.hpp"
#include "tmp/vector.hpp"

#include <type_traits>
#include <iostream>

template <class T>
struct call_tester
{
	struct type
	{
		void operator () (int const &, int const &)
		{
			std::cout << T::value << ' ';
		}
	};
};

template <class T, class U>
struct d_call_tester
{
	struct type
	{
		void operator () (int const &, int const &)
		{
			std::cout << "(" << T::value << ", " << U::value << ") ";
		}
	};
};

template <class T>
struct until_tester
{
	struct type
	{
		bool operator () (int const &)
		{
			std::cout << T::value << ' ';
			return T::value == 0;
		}
	};
};


int main(void)
{
	typedef tmp::vector<
		std::integral_constant<int, 1>,
		std::integral_constant<int, 2>
	> v1;
	
	typedef tmp::vector<
		std::integral_constant<int, -3>,
		std::integral_constant<int, 6>,
		std::integral_constant<int, 0>,
		std::integral_constant<int, -5>
	> v2;
	

	typedef call_tester<tmp::_1> call_test_t;
	typedef d_call_tester<tmp::_1, tmp::_2> d_call_test_t;
	typedef until_tester<tmp::_1> until_test_t;
	
	// Print test vectors
	std::cout << "Test vectors" << std::endl;
	std::cout << "============" << std::endl;
	std::cout << "v1 = [1, 2]" << std::endl;
	std::cout << "v2 = [-3 6 0 -5]" << std::endl;
	std::cout << std::endl;
	
	// Test call_each
	std::cout << "Testing call_each" << std::endl;
	std::cout << "=================" << std::endl;
	
	std::cout << "v1 = ";
	tmp::call_each<v1, call_test_t, int const &, int const &>(0, 0);
	std::cout << std::endl;
	
	std::cout << "v2 = ";
	tmp::call_each<v2, call_test_t, int const &, int const &>(0, 0);
	std::cout << std::endl;
	std::cout << std::endl;
	
	// Test d_call_each
	std::cout << "Testing d_call_each" << std::endl;
	std::cout << "===================" << std::endl;
	
	std::cout << "v1 x v2 = ";
	tmp::d_call_each<v1, v2, d_call_test_t, int const &, int const&>(0, 0);
	std::cout << std::endl;
	
	std::cout << "v2 x v1 = ";
	tmp::d_call_each<v2, v1, d_call_test_t, int const &, int const&>(0, 0);
	std::cout << std::endl;
	std::cout << std::endl;
	
	// Test call_until
	std::cout << std::boolalpha;
	std::cout << "Testing call_until" << std::endl;
	std::cout << "==================" << std::endl;
	std::cout << "v1 == 0? : ";
	std::cout << tmp::call_until<v1, until_test_t, int const&>(0) << std::endl;
	std::cout << "v2 == 0? : ";
	std::cout << tmp::call_until<v2, until_test_t, int const&>(0) << std::endl;
	std::cout << std::endl; 
	
	return 0;
}
