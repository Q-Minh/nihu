/**
* \file element.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief Declaration of class ::element and its specialisations
*/
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include "shapeset.hpp"

/**
 * \brief class describing the overlapping state of two element
 */
class element_overlapping
{
	unsigned num;			/**< \brief number of coincident nodes */
	unsigned ind1, ind2;	/**< \brief start node indices */
public:
	/** \brief return number of coincident nodes */
	unsigned get_num(void) const  {	return num; }
	/** \brief return index of first coincident node in element 1 */
	unsigned get_ind1(void) const { return ind1; }
	/** \brief return index of first coincident node in element 2 */
	unsigned get_ind2(void) const { return ind2; }

	/** \brief constructor */
	element_overlapping(unsigned num = 0, unsigned ind1 = 0, unsigned ind2 = 0) :
		num(num), ind1(ind1), ind2(ind2)
	{
	}
};


template <class Derived>
struct element_traits;

template <class Derived>
class element_base
{
public:
	typedef element_traits<Derived> traits_t;

	/** \brief the dimension of the element's location variable \f$x\f$ */
	static unsigned const x_dim = traits_t::x_dim;
	/** \brief the elements's L-set */
	typedef typename traits_t::lset_t lset_t;

	/** \brief the elements's domain */
	typedef typename lset_t::domain_t domain_t;
	/** \brief the base domain's scalar variable */
	typedef typename domain_t::scalar_t scalar_t;
	/** \brief the dimension of the element's domain variable \f$\xi\f$ */
	static unsigned const xi_dim = domain_t::dimension;

	/** \brief number of shape functions in the set, inherited from the LSet */
	static unsigned const num_nodes = lset_t::num_nodes;
	/** \brief type of the shape functions' independent variable \f$\xi\f$ */
	typedef typename lset_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector, inherited from the LSet */
	typedef typename lset_t::shape_t L_t;
	/** \brief type of an \f$\nabla L(\xi)\f$ gradient matrix */
	typedef typename lset_t::dshape_t dL_t;

	/** \brief type of the element's independent location variable \f$x\f$ */
	typedef Eigen::Matrix<scalar_t, 1, x_dim> x_t;
	/** \brief type of the gradient of the element's independent location variable \f$x'_{\xi}\f$ */
	typedef Eigen::Matrix<scalar_t, xi_dim, x_dim> dx_t;
	/** \brief matrix type that stores the element's corner coordinates \f$x_i\f$ */
	typedef Eigen::Matrix<unsigned, 1, num_nodes> nodes_t;
	/** \brief matrix type that stores the element's corner coordinates \f$x_i\f$ */
	typedef Eigen::Matrix<scalar_t, num_nodes, x_dim> coords_t;
	/** \brief type that stores the element's id */
	typedef Eigen::Matrix<unsigned, 1, 1> id_t;


protected:
	/** \brief the element's identifier */
	id_t m_id;
	/** \brief the element's nodal indices in the mesh */
	nodes_t m_nodes;
	/** \brief the element's corner coordinates \f$x_i\f$ */
	coords_t m_coords;
	/** \brief the element's ceter */
	x_t m_center;

public:
	/**
	* \brief default constructor for std::vector container
	*/
	element_base() {}

	/**
	* \brief constructor
	* \param [in] id the elem id
	* \param [in] nodes the elem node indices
	* \param [in] coords the elem coordinates
	*/
	element_base(unsigned id, nodes_t const &nodes, coords_t const &coords)
		: m_nodes(nodes), m_coords(coords), m_center(get_x(domain_t::get_center()))
	{
		this->m_id << id;
	}

	/**
	* \brief return elem id
	* \return elem id
	*/
	id_t const &get_id(void) const
	{
		return m_id;
	}

	/**
	* \brief return nodes
	* \return elem nodes
	*/
	nodes_t const &get_nodes(void) const
	{
		return m_nodes;
	}

	/**
	* \brief return coordinates
	* \return elem coordinates
	*/
	coords_t const &get_coords(void) const
	{
		return m_coords;
	}

	/**
	* \brief return element location
	* \param [in] \xi location \f$\xi\f$ in the base domain
	* \return location \f$x\f$ in the element
	*/
	x_t get_x(xi_t const &xi) const
	{
		return lset_t::eval_shape(xi).transpose() * m_coords;
	}

	/**
	* \brief return element center
	* \return location \f$x_0\f$ in the element center
	*/
	x_t const &get_center(void) const
	{
		return m_center;
	}

	/**
	* \brief return element location gradient
	* \param [in] \xi location \f$\xi\f$ in the base domain
	* \return location gradient \f$x'_{\xi}\f$ in the element
	*/
	dx_t get_dx(xi_t const &xi) const
	{
		return lset_t::eval_dshape(xi).transpose() * m_coords;
	}

	/**
	 * \brief check overlapping state with other element
	 * \tparam OtherElem type of the other element
	 * \param [in] other reference to the other element
	 * \return structure containing overlapping information
	 */
	template <class OtherElem>
	element_overlapping get_overlapping(OtherElem const &other) const
	{
		unsigned num_coinc = 0;	// number of coincident nodes
		unsigned start_ind1 = 0, start_ind2 = 0;

		auto const &otherNodes = other.get_nodes();
		for (unsigned i = 0; i < num_nodes; ++i)
		{
			for (unsigned j = 0; j < OtherElem::num_nodes; ++j)
			{
				if (m_nodes[i] == otherNodes[j])
				{
					// no previous coincident nodes
					if (num_coinc == 0)
					{
						start_ind1 = i;
						start_ind2 = j;
					}
					else
					{
						if (i < start_ind1
							|| (i==num_nodes-1
									&& start_ind1 == 0
									&& num_coinc < num_nodes-1))
							start_ind1 = i;
						if (j < start_ind2 ||
							(j==OtherElem::num_nodes-1 &&
								start_ind2 == 0 &&
								num_coinc < OtherElem::num_nodes-1))
							start_ind2 = j;
					}
					num_coinc++;
				}
			}
		}
		if (num_coinc == 0)
			return element_overlapping();
		return element_overlapping(num_coinc, start_ind1, start_ind2);
	}
};


class tria_1_elem;

template <>
struct element_traits<tria_1_elem>
{
	static const unsigned x_dim = 3;
	typedef tria_1_shape_set lset_t;
};

class tria_1_elem : public element_base<tria_1_elem>
{
public:
	tria_1_elem(unsigned id, nodes_t const &nodes, coords_t const &coords)
		: element_base(id, nodes, coords)
	{
		m_normal = (coords.row(1)-coords.row(0)).cross(coords.row(2)-coords.row(0));
	}

	x_t const &get_normal(xi_t const &) const
	{
		return m_normal;
	}

protected:
	x_t m_normal;
};


class quad_1_elem;

template <>
struct element_traits<quad_1_elem>
{
	static const unsigned x_dim = 3;
	typedef quad_1_shape_set lset_t;
};

class quad_1_elem : public element_base<quad_1_elem>
{
public:
	quad_1_elem(unsigned id, nodes_t const &nodes, coords_t const &coords)
		: element_base(id, nodes, coords)
	{
		x_t a = (coords.row(1)-coords.row(0)+coords.row(2)-coords.row(3))/4.0;
		x_t b = (coords.row(0)-coords.row(1)+coords.row(2)-coords.row(3))/4.0;
		x_t c = (coords.row(2)-coords.row(1)+coords.row(3)-coords.row(0))/4.0;

		m_n0 = a.cross(c);
		m_n_xi.row(0) = a.cross(b);
		m_n_xi.row(1) = b.cross(c);
	}

	x_t get_normal(xi_t const &xi) const
	{
		return m_n0 + xi.transpose() * m_n_xi;
	}

protected:
	x_t m_n0;
	Eigen::Matrix<double, 2, 3> m_n_xi;
};

#endif

