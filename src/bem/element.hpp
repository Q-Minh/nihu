/**
* \file element.hpp
* \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
* \brief Declaration of element classes and specialisations
*/
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include "space.hpp"
#include "shapeset.hpp"

/**
* \brief class describing the overlapping state of two element
*/
class element_overlapping
{
public:
	/** \brief return number of coincident nodes */
	unsigned get_num(void) const
	{
		return num;
	}
	
	/** \brief return index of first coincident node in element 1 */
	unsigned get_ind1(void) const
	{
		return ind1;
	}
	
	/** \brief return index of first coincident node in element 2 */
	unsigned get_ind2(void) const
	{
		return ind2;
	}

	/** \brief constructor */
	element_overlapping(unsigned num = 0, unsigned ind1 = 0, unsigned ind2 = 0) :
		num(num), ind1(ind1), ind2(ind2)
	{
	}

private:
	/** \brief number of coincident nodes */
	unsigned num;
	/** \brief start node indices */
	unsigned ind1, ind2;
};


/** \brief Traits class describing element properties */
template <class Derived>
struct element_traits;

/** \brief Metafunction providing an element id */
template <class elem_t>
struct elem_id
{
	static unsigned const value =
		element_traits<elem_t>::space_t::dimension * 10000 +
		shape_set_id<typename element_traits<elem_t>::lset_t>::value;
};

/**
* \brief The geometrical element representation
* \tparam Derived CRTP derived class
*/
template <class Derived>
class element_base
{
public:
	/** \brief the CRTP derived class */
	typedef element_traits<Derived> traits_t;

	/** \brief the element space type */
	typedef typename traits_t::space_t space_t;
	/** \brief the dimension of the element's location variable \f$x\f$ */
	static unsigned const x_dim = space_t::dimension;

	/** \brief the elements's L-set */
	typedef typename traits_t::lset_t lset_t;

	/** \brief the element id */
	static unsigned const id = elem_id<Derived>::value;

	/** \brief the elements's reference domain */
	typedef typename lset_t::domain_t domain_t;
	/** \brief the base domain's scalar variable */
	typedef typename space_t::scalar_t scalar_t;
	/** \brief the dimension of the element's domain variable \f$\xi\f$ */
	static unsigned const xi_dim = domain_t::dimension;

	/** \brief number of shape functions in the L-set */
	static unsigned const num_nodes = lset_t::num_nodes;
	/** \brief type of the shape functions' independent variable \f$\xi\f$ */
	typedef typename lset_t::xi_t xi_t;
	/** \brief type of an \f$L(\xi)\f$ vector */
	typedef typename lset_t::shape_t L_t;
	/** \brief type of an \f$\nabla L(\xi)\f$ gradient matrix */
	typedef typename lset_t::dshape_t dL_t;

	/** \brief type of the element's independent location variable \f$x\f$ */
	typedef typename space_t::location_t x_t;
	/** \brief type of the gradient of the element's independent location variable \f$x'_{\xi}\f$ */
	typedef Eigen::Matrix<scalar_t, x_dim, xi_dim> dx_t;
	/** \brief matrix type that stores the element's corner nodes \f$x_i\f$ */
	typedef Eigen::Matrix<unsigned, num_nodes, 1> nodes_t;
	/** \brief matrix type that stores the element's corner coordinates \f$x_i\f$ */
	typedef Eigen::Matrix<scalar_t, x_dim, num_nodes> coords_t;
	/** \brief type that stores the element's id */
	typedef Eigen::Matrix<unsigned, 1, 1> id_t;


protected:
	/** \brief the element's identifier */
	id_t m_id;
	/** \brief the element's nodal indices in the mesh */
	nodes_t m_nodes;
	/** \brief the element's corner coordinates \f$x_i\f$ */
	coords_t m_coords;
	/** \brief the element's center */
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
	element_base(coords_t const &coords, unsigned id = 0, nodes_t const &nodes = nodes_t())
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
		return m_coords * lset_t::eval_shape(xi);
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
		return m_coords * lset_t::eval_dshape(xi);
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


/** \brief a linear triangle element in 3D space */
class tria_1_elem;

/** \brief element properties of a linear 3D tria element */
template <>
struct element_traits<tria_1_elem>
{
	typedef space_3d space_t;
	/** \brief the shape set */
	typedef tria_1_shape_set lset_t;
};

/** \brief a linear tria element in 3D space */
class tria_1_elem : public element_base<tria_1_elem>
{
public:
	/**
	* \brief constructor
	* \param [in] id the element id
	* \param [in] nodes the nodal indices
	* \param [in] coords the nodal coordinates
	*/
	tria_1_elem(coords_t const &coords, unsigned id = 0, nodes_t const &nodes = nodes_t())
		: element_base(coords, id, nodes)
	{
		dx_t dx0(get_dx(xi_t::Zero()));
		m_normal = dx0.col(0).cross(dx0.col(1));
	}

	/** \brief return normal vector
	* \return element normal
	*/
	x_t const &get_normal(xi_t const &) const
	{
		return m_normal;
	}

protected:
	/** \brief the normal vector */
	x_t m_normal;
};


/** \brief a linear quad element in 3D space */
class quad_1_elem;

/** \brief element properties of a linear 3D quad element */
template <>
struct element_traits<quad_1_elem>
{
	/** \brief dimensionality of the x space */
	typedef space_3d space_t;
	/** \brief the shape set */
	typedef quad_1_shape_set lset_t;
};

/** \brief a linear quad element in 3D space */
class quad_1_elem : public element_base<quad_1_elem>
{
public:
	/**
	* \brief constructor
	* \param [in] id the element id
	* \param [in] nodes the nodal indices
	* \param [in] coords the nodal coordinates
	*/
	quad_1_elem(coords_t const &coords, unsigned id = 0, nodes_t const &nodes = nodes_t())
		: element_base(coords, id, nodes)
	{
		dx_t dx0 = get_dx(xi_t::Zero());
		dx_t dx1 = get_dx(xi_t::Ones())-dx0;

		m_n0 = dx0.col(0).cross(dx0.col(1));
		m_n_xi.col(0) = dx0.col(0).cross(dx1.col(0));
		m_n_xi.col(1) = dx1.col(0).cross(dx0.col(1));
	}

	/** \brief return normal vector
	* \param [in] xi local coordinate
	* \return element normal
	*/
	x_t get_normal(xi_t const &xi) const
	{
		return m_n0 + m_n_xi * xi;
	}

protected:
	/** \brief the normal at xi=0 */
	x_t m_n0;
	/** \brief 1st order tag of the Taylor series of the normal */
	dx_t m_n_xi;
};



// forward declaration
template <class LSet, class scalar_t>
class general_surface_element;

/** \brief specialisation of element_traits for the general surface element */
template <class LSet, class scalar_t>
struct element_traits<general_surface_element<LSet, scalar_t> >
{
	/** \brief dimensionality of the x space */
	typedef space<scalar_t, LSet::domain_t::dimension + 1> space_t;
	/** \brief the shape set */
	typedef LSet lset_t;
};

/** \brief class describing a general surface element that computes its normal in the general way
* \tparam LSet type of the geometry shape set
*/
template <class LSet, class scalar_t>
class general_surface_element : public element_base<general_surface_element<LSet, scalar_t> >
{
public:
	/** \brief the base class type */
	typedef element_base<general_surface_element<LSet, scalar_t> > base_t;

	/** \brief type of the coordinate matrix */
	typedef typename base_t::coords_t coords_t;
	/** \brief type of the node index vector */
	typedef typename base_t::nodes_t nodes_t;
	/** \brief type of the coordinate  on the standard element */
	typedef typename base_t::xi_t xi_t;
	/** \brief type of the coordinate vector */
	typedef typename base_t::x_t x_t;
	/** \brief type of the coordinate derivative vector */
	typedef typename base_t::dx_t dx_t;

	using base_t::get_dx;

	/** \brief constructor
	* \param [in] coords the coordinate matrix
	* \param [in] id element id
	* \param [in] nodes the nodal index vector
	*/
	general_surface_element(coords_t const &coords, unsigned id = 0, nodes_t const &nodes = nodes_t())
		: base_t(coords, id, nodes)
	{
	}

	/** \brief return normal vector at given location
	* \param [in] xi the location in the standard domain
	* \return normal vector
	*/
	x_t get_normal(xi_t const &xi) const
	{
		dx_t dx = get_dx(xi);
		return dx.col(0).cross(dx.col(1));
	}
};

/** \brief quadratic 6-noded triangle element */
typedef general_surface_element<tria_2_shape_set, tria_1_elem::space_t::scalar_t> tria_2_elem;
/** \brief quadratic 9-noded quadrilateral element */
typedef general_surface_element<quad_2_shape_set, quad_1_elem::space_t::scalar_t> quad_2_elem;

#endif

