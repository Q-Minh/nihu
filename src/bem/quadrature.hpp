/**
 * \file quadrature.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief implementation of class ::quadrature_elem, ::quadrature_base
 */
#ifndef QUADRATURE_HPP_INCLUDED
#define QUADRATURE_HPP_INCLUDED

#include "shapeset.hpp"

#include <type_traits> // is_same

#include <iostream>
#include <stdexcept>

#include <numeric>

#include <Eigen/StdVector>
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >


/**
 * \brief a quadrature element is a base point and a weight
 * \tparam XiType the local coordinate type
 * \tparam ScalarType the scalar of given coordinates
 */
template <class XiType, class ScalarType>
class quadrature_elem
{
public:
	typedef XiType xi_t;			/**< \brief template argument as tested type */
	typedef ScalarType scalar_t;	/**< \brief template argument as nested type */

	/**
	 * \brief constructor initialising all members
	 * \param [in] xi base location
	 * \param [in] w weight
	 */
	quadrature_elem(xi_t const &xi = xi_t(), scalar_t const &w = scalar_t())
		: m_xi(xi), m_w(w)
	{
	}

	/**
	 * \brief return constant reference to base point
	 * \return reference to base point
	 */
	xi_t const &get_xi(void) const
	{
		return m_xi;
	}

	/**
	 * \brief return constant reference to weight
	 * \return reference to weight
	 */
	scalar_t const &get_w(void) const
	{
		return m_w;
	}

	/**
	 * \brief transform a location according to a shape function set
	 * \tparam LSet the shape function set
	 * \param [in] coords transformation corners
	 * \return reference to the transformed object
	 */
	template <class LSet>
	quadrature_elem &transform_inplace(const Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> &coords)
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// compute Jacobian and multiply with original weight
		m_w *= (LSet::eval_dshape(m_xi).transpose() * coords).determinant();
		// compute new quadrature location
		m_xi = LSet::eval_shape(m_xi).transpose() * coords;
		return *this;
	}

	// struct is used in a std::vector, therefore this declaration is necessary
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	xi_t m_xi;		/**< \brief base point */
	scalar_t m_w;	/**< \brief weight */
};


// forward declaration of singular traits
template <class quadrature>
struct singular_traits;

// forward declaration of quadrature traits
template <class Derived>
struct quadrature_traits;

/**
 * \brief CRTP base class of all quadratures
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
class quadrature_base :
	public std::vector<quadrature_elem<typename quadrature_traits<Derived>::domain_t::xi_t, typename quadrature_traits<Derived>::domain_t::scalar_t> , Eigen::aligned_allocator<quadrature_elem<typename quadrature_traits<Derived>::domain_t::xi_t, typename quadrature_traits<Derived>::domain_t::scalar_t> > >
{
protected:
	/**
	 * \brief static cast to derived const type
	 */
	Derived const &derived(void) const
	{
		return static_cast<Derived const &>(*this);
	}

	/**
	 * \brief static cast to derived type
	 */
	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

public:
	typedef quadrature_traits<Derived> traits_t;				/**< \brief corresponding traits type */
	typedef typename traits_t::domain_t domain_t;				/**< \brief the domain type */
	typedef typename domain_t::xi_t xi_t; 						/**< \brief local coordinate type */
	typedef typename domain_t::scalar_t scalar_t; 				/**< \brief local scalar type */
	typedef quadrature_elem<xi_t, scalar_t> quadrature_elem_t;	/**< \brief the quadrature elem type */

	/**
	 * \brief constructor allocating space for the quadrature elements
	 * \param N number of quadrature elements
	 */
	quadrature_base(unsigned N = 0)
	{
		this->reserve(N);
	}

	/**
	 * \brief print a quadrature into an output stream
	 * \param os the output stream
	 * \return reference to the stream
	 */
	std::ostream & print (std::ostream & os) const
	{
		os << "Base points: \t Weights:" << std::endl;
		std::for_each(this->begin(), this->end(), [&os] (quadrature_elem_t const &e) {
			os << e.get_xi().transpose() << "\t\t" << e.get_w() << std::endl;
		});
		os << "Sum of weights: " <<
			std::accumulate(this->begin(), this->end(), scalar_t(), [] (scalar_t x, quadrature_elem_t const &qe) {
			return x + qe.get_w();
		}) << std::endl;

		return os;
	}

	/**
	 * \brief transform the domain of the quadrature with a given shape set and corner points
	 * \tparam LSet shape set type of transformation
	 * \param [in] coords corner coordinates of transformed domain
	 * \return transformed quadrature
	 */
	template <class LSet>
	Derived transform(const Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> &coords) const
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// dimension check
		static_assert(std::is_same<xi_t, typename LSet::xi_t>::value,
			"Quadrature and shape set dimensions must match");
		Derived result(derived());
		return result.transform_inplace<LSet>(coords);
	}

	/**
	 * \brief transform the domain of the quadrature in place
	 * \tparam LSet shape set type of transformation
	 * \param [in] coords corner coordinates of transformed domain
	 * \return reference to the transformed quadrature
	 */
	template <class LSet>
	Derived &transform_inplace(const Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> &coords)
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// dimension check
		static_assert(std::is_same<xi_t, typename LSet::xi_t>::value,
			"Quadrature and shape set dimensions must match");
		for (auto it = this->begin(); it != this->end(); ++it)
			it->transform_inplace<LSet>(coords);
		return derived();
	}

	/**
	 * \brief add two quadratures
	 * \details It is assured that the dimensions match.
	 * \tparam otherDerived other quadrature type
	 * \param [in] other the other quadrature
	 * \return new quadrature
	 */
	template <class otherDerived>
	Derived operator +(const quadrature_base<otherDerived> &other) const
	{
		static_assert(std::is_same<xi_t, typename otherDerived::xi_t>::value,
			"Quadrature domain dimensions must match");
		Derived result(this->size()+other.size());
		result.insert(result.end(), this->begin(), this->end());
		result.insert(result.end(), other.begin(), other.end());
		return result;
	}

	/**
	 * \brief add another quadrature to this
	 * \details It is assured that the dimensions match.
	 * \tparam otherDerived other quadrature type
	 * \param [in] other the other quadrature
	 * \return reference to the new quadrature
	 */
	template <class otherDerived>
	Derived &operator +=(const quadrature_base<otherDerived> &other)
	{
		static_assert(std::is_same<xi_t, typename otherDerived::xi_t>::value,
			"Quadrature domain dimensions must match");
		this->reserve(this->size()+other.size());
		this->insert(this->end(), other.begin(), other.end());
		return derived();
	}

	/**
	 * \brief create a singular quadrature from a selected internal point
	 * \param [in] degree quadrature degree of the generating element
	 * \param [in] singular_point the internal singular point
	 */
	static Derived singular_quadrature_inside(unsigned degree, xi_t const &singular_point)
	{
		typedef singular_traits<Derived> singular_traits_t;
		typedef typename singular_traits_t::singular_source_type singular_source_t;
		typedef typename singular_traits_t::transformation_lset transformation_lset;
		const unsigned nResolution = singular_traits_t::nResolution;
		const unsigned nCorners = singular_traits_t::nCorners;

		xi_t const *corners = domain_t::get_corners();

		Derived result;
		singular_source_t source_quad(degree);
		for (unsigned r = 0; r < nResolution; ++r)
		{
			Eigen::Matrix<scalar_t, nCorners, domain_t::dimension> coords;
			for (unsigned c = 0; c < nCorners; ++c)
			{
				if (singular_traits_t::index_inside[r][c] == -1)
					coords.row(c) = singular_point;
				else
					coords.row(c) = corners[singular_traits_t::index_inside[r][c]];
			}
			result += source_quad.transform<transformation_lset>(coords);
		}
		return result;
	}
}; // quadrature_base


/**
 * \brief print a quadrature into an ouput stream
 * \tparam the quadrature's domain type
 * \param os the output stream
 * \param Q the quadrature
 * \return the modified output stream
 */
template<class Derived>
std::ostream & operator << (std::ostream & os, const quadrature_base<Derived>& Q)
{
	return Q.print(os);
}

/**
 * \brief metafunction to assign a quadrature type to a quadrature family and a domain
 */
template <class Family, class Domain>
struct quadrature_type;

#endif

