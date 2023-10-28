// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/// \file quadrature.hpp
/// \ingroup quadrature
/// \brief implementation of class NiHu::quadrature_elem, NiHu::quadrature_base
#ifndef QUADRATURE_HPP_INCLUDED
#define QUADRATURE_HPP_INCLUDED

#include "shapeset.hpp"
#include "../library/lib_shape.hpp"
#include "../util/eigen_utils.hpp"
#include "../util/crtp_base.hpp"

#include <Eigen/Core>

namespace NiHu
{

/**
* \brief a quadrature element is a base point and a weight
* \tparam XiType the local coordinate type
* \tparam ScalarType the scalar of coordinates
*/
template <class XiType, class ScalarType>
class quadrature_elem
{
public:
	/** \brief template argument as tested type */
	typedef XiType xi_t;
	/** \brief template argument as nested type */
	typedef ScalarType scalar_t;

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
	 * \brief return reference to base point
	 */
	xi_t &get_xi(void)
	{
		return m_xi;
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
	* \brief set quadrature weight
	* \param [in] w new weight
	*/
	void set_w(scalar_t const &w)
	{
		m_w = w;
	}

	/**
	* \brief multiply a quadrature element by a scalar
	* \param c the scalar multiplier
	*/
	quadrature_elem &operator *=(scalar_t const &c)
	{
		m_w *= c;
		return *this;
	}

	/**
	* \brief transform a quadrature element according to a shape function set and corner nodes
	* \tparam LSet the shape function set
	* \tparam Signed indicates if the sign of the Jacobian needs to be taken into ancount
	* \param [in] coords transformation corners
	* \return reference to the transformed object
	*/
	template <class LSet, bool Signed = true>
	quadrature_elem &transform_inplace(
		Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> const &coords)
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// compute Jacobian and multiply with original weight
		auto j = (LSet::template eval_shape<1>(m_xi).transpose() * coords).determinant();
		if (Signed)
			m_w *= j;
		else
			m_w *= std::abs(j);
		// compute new quadrature location
		m_xi = LSet::template eval_shape<0>(m_xi).transpose() * coords;
		return *this;
	}

	// struct is used in a std::vector, therefore this declaration is necessary
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	/** \brief base point */
	xi_t m_xi;
	/** \brief weight */
	scalar_t m_w;
};


// forward declaration of quadrature traits
template <class Derived>
struct quadrature_traits;

/** \brief metafunction to assign a quadrature element to a quadrature */
template <class Derived>
struct quadr_elem
{
	typedef quadrature_elem<
		typename quadrature_traits<Derived>::domain_t::xi_t,
		typename quadrature_traits<Derived>::domain_t::scalar_t
	> type;
};

/**
* \brief CRTP base class of all quadratures
* \tparam Derived the CRTP derived class
*/
template <class Derived>
class quadrature_base :
	public eigen_std_vector<typename quadr_elem<Derived>::type>::type
{
public:
	/** \brief the base vector class of the quadrature */
	typedef typename eigen_std_vector<typename quadr_elem<Derived>::type>::type base_t;

	NIHU_CRTP_HELPERS

public:
	/** \brief traits type */
	typedef quadrature_traits<Derived> traits_t;
	/** \brief domain type */
	typedef typename traits_t::domain_t domain_t;
	/** \brief local coordinate type */
	typedef typename domain_t::xi_t xi_t;
	/** \brief local scalar type */
	typedef typename domain_t::scalar_t scalar_t;
	/** \brief quadrature elem type */
	typedef quadrature_elem<xi_t, scalar_t> quadrature_elem_t;

	/**
	* \brief constructor allocating space for the quadrature elements
	* \param N number of quadrature elements
	*/
	quadrature_base(size_t N = 0)
	{
		base_t::reserve(N);
	}

	/**
	* \brief return sum of quadrature weights
	* \return sum of weights
	*/
	scalar_t sum_of_weights(void) const
	{
		scalar_t res = 0.0;
		for (auto it = base_t::begin(); it != base_t::end(); ++it)
			res += it->get_w();
		return res;
	}

	/**
	* \brief print a quadrature into an output stream
	* \param os the output stream
	* \return reference to the stream
	*/
	std::ostream & print(std::ostream & os) const
	{
		os << base_t::size() << std::endl;
		for (auto it = base_t::begin(); it != base_t::end(); ++it)
			os << it->get_xi().transpose() << "\t\t" << it->get_w() << std::endl;
		return os;
	}

	/**
	* \brief multiply the quadrature by a scalar
	* \param [in] c the scalar multiplier
	* \return reference to the the new quadrature
	*/
	Derived &operator *=(scalar_t const &c)
	{
		for (size_t i = 0; i < base_t::size(); ++i)
			(*this)[i] *= c;
		return derived();
	}

	/**
	* \brief transform the domain of the quadrature with a given shape set and corner points
	* \tparam LSet shape set type of transformation
	* \tparam Signed indicates if the sign of the Jacobian needs to be taken into ancount
	* \param [in] coords corner coordinates of transformed domain
	* \return transformed quadrature
	* \todo this transformation is sick, does not fit into our quadrature scheme conceptually
	*/
	template <class LSet, bool Signed = true>
	Derived transform(
		Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> const &coords) const
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// dimension check
		static_assert(std::is_same<xi_t, typename LSet::xi_t>::value,
			"Quadrature and shape set dimensions must match");
		Derived result(derived());
		return result.template transform_inplace<LSet, Signed>(coords);
	}

	/**
	* \brief transform the domain of the quadrature in place
	* \tparam LSet shape set type of transformation
	* \tparam Signed indicates if the sign of the Jacobian needs to be taken into ancount
	* \param [in] coords corner coordinates of transformed domain
	* \return reference to the transformed quadrature
	*/
	template <class LSet, bool Signed = true>
	Derived &transform_inplace(
		const Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> &coords)
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// dimension check
		static_assert(std::is_same<xi_t, typename LSet::xi_t>::value,
			"Quadrature and shape set dimensions must match");
		for (auto it = base_t::begin(); it != base_t::end(); ++it)
			it->template transform_inplace<LSet, Signed>(coords);
		return derived();
	}

	/**
	* \brief add two quadratures
	* \tparam otherDerived other quadrature type
	* \param [in] other the other quadrature
	* \return new quadrature
	* \details It is assured that the dimensions match.
	*/
	template <class otherDerived>
	Derived operator +(const quadrature_base<otherDerived> &other) const
	{
		static_assert(std::is_same<xi_t, typename otherDerived::xi_t>::value,
			"Quadrature domain dimensions must match");
		Derived result(base_t::size()+other.size());
		result.insert(result.end(), base_t::begin(), base_t::end());
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
		base_t::reserve(base_t::size()+other.size());
		base_t::insert(base_t::end(), other.begin(), other.end());
		return derived();
	}
}; // class quadrature_base

/**
* \brief print a quadrature into an ouput stream
* \tparam Derived the quadrature's type
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

}


#endif
