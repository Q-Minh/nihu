#ifndef NIHU_JACOBIAN_COMPUTER_HPP_INCLUDED
#define NIHU_JACOBIAN_COMPUTER_HPP_INCLUDED

namespace NiHu
{

template <class elem, class enable = void>
class jacobian_computer;

template <class elem_t>
class jacobian_computer<
	elem_t,
	typename std::enable_if<
	NiHu::element_traits::is_surface_element<elem_t>::value
	>::type
>
{
public:
	static double eval(elem_t const &elem, typename elem_t::xi_t const &xi)
	{
		return elem.get_normal(xi).norm();
	}
};

template <class elem_t>
class jacobian_computer<
	elem_t,
	typename std::enable_if<
	!NiHu::element_traits::is_surface_element<elem_t>::value
	>::type
>
{
public:
	static double eval(elem_t const &elem, typename elem_t::xi_t const &xi)
	{
		return elem.get_dx(xi).determinant();
	}
};

}// namespace NiHu

#endif // JACOBIAN_COMPUTER_INCLUDED
