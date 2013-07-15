
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
			std::cout << "location ctor\n";
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
struct normal_jacobian
{
	template <class wall>
	struct brick : public wall
	{
	public:
		typedef Space space_t;
		typedef typename space_t::location_t x_t;
		typedef typename space_t::scalar_t scalar_t;

		template <class elem_t>
		brick(elem_t const &elem, typename elem_t::xi_t const &xi) :
			wall(elem, xi),
			m_norm(elem.get_normal(xi)),
			m_jac(m_norm.norm())
		{
			m_norm /= m_jac;
			std::cout << "normal ctor\n";
		}

		x_t const &get_normal(void) const
		{
			return m_norm;
		}

		scalar_t const &get_jacobian(void) const
		{
			return m_jac;
		}

	private:
		x_t m_norm;
		scalar_t m_jac;
	};
};

#endif // LOCATION_HPP_INCLUDED

