#include "tmp/defines_type.hpp"
#include <iostream>

// metafunction declarations with macro
DEFINES_TYPE(default_version)
DEFINES_TYPE(special_version)

template <class C>
struct test
{
	typedef void default_version;
};

template <>
struct test<bool>
{
	typedef bool special_version;
};

int main(void)
{
	std::cout << std::boolalpha;

	std::cout << defines_default_version<test<void> >::value << ' '
		<< defines_special_version<test<void> >::value << std::endl;

	std::cout << defines_default_version<test<bool> >::value << ' '
		<< defines_special_version<test<bool> >::value << std::endl;
}
