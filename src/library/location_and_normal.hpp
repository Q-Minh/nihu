/**
 * \file location_and_normal.hpp
 * \ingroup library
 * \brief implementation of kernel inputs ::location and ::location_with_normal
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */

#ifndef LOCATION_AND_NORMAL_HPP_INCLUDED
#define LOCATION_AND_NORMAL_HPP_INCLUDED

#include "../bem/kernel_input.hpp"

// forward declaration
template <class space_t>
class location;

/** \brief traits of a location
* \tparam Space the x-space
*/
template <class Space>
struct kernel_input_traits<location<Space> >
{
	/** \brief the space type */
	typedef Space space_t;
};


/** \brief a kernel input representing a single location \f$ \gamma = {\bf x} \f$
* \tparam Space the coordinate space of the location
*/
template <class Space>
class location :
	public kernel_input_base<location<Space> >
{
public:
	/** \brief the CRTP base class */
	typedef kernel_input_base<location<Space> > base_t;
	/** \brief the space type */
	typedef typename base_t::space_t space_t;
	/** \brief the location vector type */
	typedef typename space_t::location_t x_t;

	/** \brief constructor from element and reference domain vector
	* \param [in] elem the element to construct from
	* \param [in] the location in the reference domain
	*/
	template <class elem_t>
	location(elem_t const &elem, typename elem_t::xi_t const &xi)
		: base_t(elem, xi), m_x(elem.get_x(xi))
	{
	}

	x_t const &get_x(void) const
	{
		return m_x;
	}

protected:
	x_t m_x;	/**< \brief the location */
};


// forward declaration
template <class space_t>
class location_with_normal;


/** \brief traits of a location with normal
* \tparam Space the x-space
*/
template <class Space>
struct kernel_input_traits<location_with_normal<Space> >
{
	typedef Space space_t;
};


/** \brief a kernel input representing a location and a normal \f$ \gamma = \left\{{\bf x}, {\bf n}({\bf x})\right\} \f$
* \tparam Space the coordinate space of the location
*/
template <class Space>
class location_with_normal :
	public kernel_input_base<location_with_normal<Space> >,
	public jacobian<typename Space::scalar_t>
{
public:
	typedef kernel_input_base<location_with_normal<Space> > base_t;
	typedef typename base_t::space_t space_t;
	typedef typename space_t::location_t x_t;

	template <class elem_t>
	location_with_normal(elem_t const &elem, typename elem_t::xi_t const &xi)
		: base_t(elem, xi), m_x(elem.get_x(xi))
	{
		m_unit_normal = elem.get_normal(xi);
		this->set_jacobian(m_unit_normal.norm());
		m_unit_normal /= this->get_jacobian();
	}

	x_t const &get_x(void) const
	{
		return m_x;
	}

	x_t const &get_unit_normal(void) const
	{
		return m_unit_normal;
	}

protected:
	/** \brief the location */
	x_t m_x;
	/** \brief the unit normal vector */
	x_t m_unit_normal;
};

#endif // LOCATION_AND_NORMAL_HPP_INCLUDED

