#ifndef TYPE_INFO_HPP_INCLUDED
#define TYPE_INFO_HPP_INCLUDED

#include <type_traits>
#include <iostream>

template <class T>
std::ostream & print_type_info(std::ostream &os = std::cout)
{
	os << (std::is_const<typename std::remove_reference<T>::type >::value ? "const " : "");
	os << (std::is_same<typename std::decay<T>::type, int>::value ? "int " : "");
	os << (std::is_same<typename std::decay<T>::type, char>::value ? "char " : "");
	os << (std::is_lvalue_reference<T>::value ? "&" : "");
	os << (std::is_rvalue_reference<T>::value ? "&&" : "");
	return os;
}

#endif // TYPE_INFO_HPP_INCLUDED

