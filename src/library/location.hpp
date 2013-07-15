
#ifndef LOCATION_HPP_INCLUDED
#define LOCATION_HPP_INCLUDED

#include "../bem/element.hpp"

template <class Space>
struct location
{
	template <class wall>
	class brick : public wall
	{
	public:
		typedef Space space_t;
		typedef typename space_t::location_t x_t;

		template <class elem_t>
		brick(elem_t const &elem, typename elem_t::xi_t const &xi) :
			wall(elem, xi),
			m_x(elem.get_x(xi))
		{
			std::cout << "m_x\n";
		}

		x_t const &get_x(void) const
		{
			return m_x;
		}

	private:
		x_t m_x;
	};
};


template <class Space>
struct area_elem
{
	template <class wall>
	struct brick : public wall
	{
	public:
		typedef Space space_t;
		typedef typename space_t::location_t x_t;

		template <class elem_t>
		brick(elem_t const &elem, typename elem_t::xi_t const &xi) :
			wall(elem, xi),
			m_dA(elem.get_normal(xi))
		{
			std::cout << "m_dA\n";
		}

		x_t const &get_dA(void) const
		{
			return m_dA;
		}

	private:
		x_t m_dA;
	};
};


#endif // LOCATION_HPP_INCLUDED

