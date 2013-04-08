/**
 * \file field_type_accelerator.hpp
 * \brief declaration of class field_type_accelerator
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef FIELD_TYPE_ACCELERATOR_HPP_INCLUDED
#define FIELD_TYPE_ACCELERATOR_HPP_INCLUDED

#include "field.hpp"
#include "quadrature.hpp"
#include "shapeset.hpp"

#include <vector>

/** \brief abbreviation for Eigen's std::vector declaration */
#include <Eigen/StdVector>
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >

/**
 * \brief class to store quadrature points and shape functions for a field type
 * \details field_type_accelerator is used to precompute and store quadrature points and
 * shape function values for a specified quadrature order. The class provides an iterator that
 * can traverse through the stored values
 * \tparam Field the field type that is accelerated
 */
template <class Field>
class field_type_accelerator
{
public:
	/** \brief template argument as nested type */
	typedef Field field_t;

	/** \brief the quadrature vector type */
	typedef gauss_quadrature<typename field_t::elem_t::domain_t> quadrature_t;
	/** \brief the type of one quadrature point */
	typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;
	/** \brief the quadrature iterator type */
	typedef typename quadrature_t::const_iterator quad_iter_t;

	/** \brief the shape function class's type */
	typedef typename field_t::nset_t nset_t;
	/** \brief type of a shape vector \f$N(\xi)\f$ */
	typedef typename nset_t::shape_t shape_t;
	/** \brief type of the shape function container */
	typedef EIGENSTDVECTOR(shape_t) nset_vector_t;
	/** \brief type of the shape function container's iterator */
	typedef typename nset_vector_t::const_iterator shape_iter_t;

	// forward declaration so that const_iterator can be friend of iter_tuple
	class const_iterator;

	/**
	 * \brief a pair of iterators pointing to the separate container members
	 */
	class iter_tuple
	{
		friend class const_iterator;
	public:
		/** \brief constructor from two iterators */
		iter_tuple(quad_iter_t const &qe, shape_iter_t const &sh)
			: m_quad_it(qe), m_shape_it(sh)
		{
		}

		/** \brief return the pointed quadrature element
		 * \return the pointed quadrature element
		 */
		quadrature_elem_t const &get_quadrature_elem(void) const
		{
			return *m_quad_it;
		}

		/** \brief return the pointed shape function
		 * \return the pointed shape function
		 */
		shape_t const &get_shape(void) const
		{
			return *m_shape_it;
		}

	private:
		quad_iter_t m_quad_it;	/**< \brief the quadrature iterator member */
		shape_iter_t m_shape_it;	/**< \brief the shape function iterator member */
	};

	/**
	 * \brief nested iterator class to provide access to the stored values
	 * \details instances of the class can be dereferenced by conventional iterator semantics
	 * to access the stored quadrature points and shape function values
	 */
	class const_iterator
	{
	public:
		/** \brief constructor from an iterator tuple
		 * \param [in] itpl the iterator tuple to be stored
		 */
		const_iterator (iter_tuple const &itpl) : m_tuple(itpl)
		{
		}

		/** \brief preincrement operator
		 * \param [in] the operator simply increments both stored iterators
		 */
		const_iterator &operator++()
		{
			++m_tuple.m_quad_it;
			++m_tuple.m_shape_it;
			return *this;
		}

		/** \brief inequality operator
		 * \param [in] other the other iterator to compare against
		 * \details as the two nested iterators are always incremented together, only match of
		 * the first one is checked.
		 */
		bool operator!=(const_iterator const &other)
		{
			return m_tuple.m_quad_it != other.m_tuple.m_quad_it;
		}

		/** \brief dereference operator
		 * \return an iter_tuple that can access both pointed elements
		 */
		iter_tuple const &operator*() const
		{
			return m_tuple;
		}

		/** \brief dereference operator
		 * \return address of the stored iter_tuple that can access both pointed elements
		 */
		iter_tuple const *operator->() const
		{
			return &m_tuple;
		}

	private:
		iter_tuple m_tuple;	/**< \brief a pair of iterators stored as member */
	};

	/**
	 * \brief constructor for a given quadrature order
	 * \param [in] order the quadrature order
	 * \details the constructor computes the quadrature points
	 * and the shape functions in the quadrature points
	 */
	field_type_accelerator(unsigned order) : m_quadrature(order)
	{
		m_nset_vector.reserve(m_quadrature.size());
		for (quad_iter_t quad_it = m_quadrature.cbegin(); quad_it != m_quadrature.cend(); ++quad_it)
			m_nset_vector.push_back(nset_t::eval_shape(quad_it->get_xi()));
	}

	/**
	 * \brief return begin iterator
	 * \return begin iterator of the accelerator container
	 * \details by dereferencing the begin iterator, one can access a quadrature point
	 * and a shape function as
	 * it->get_shape() or it->get_quadrature_elem
	 */
	const_iterator cbegin(void) const
	{
		return const_iterator(iter_tuple(m_quadrature.cbegin(), m_nset_vector.cbegin()));
	}

	/**
	 * \brief return end iterator
	 * \return end iterator of the accelerator container
	 */
	const_iterator cend(void) const
	{
		return const_iterator(iter_tuple(m_quadrature.cend(), m_nset_vector.cend()));
	}

protected:
	quadrature_t m_quadrature;	/**< \brief the stored quadrature points */
	nset_vector_t m_nset_vector;	/**< \brief the stored shape functions */
};


/**
 * \brief container class to store accelerators with different quadrature orders
 * \tparam Field the field type that is accelerated
 */
template <class Field>
class field_type_accelerator_pool : public std::vector<field_type_accelerator<Field> *>
{
public:
	/** \brief template argument as nested type */
	typedef Field field_t;
	/** \brief maximum order of quadratures */
	static const unsigned MAX_ORDER = 10;

	/**
	 * \brief allocate memory and initialise accelerators
	 */
	field_type_accelerator_pool(void)
	{
		this->reserve(MAX_ORDER);
		for (unsigned order = 0; order < MAX_ORDER; ++order)
			this->push_back(new field_type_accelerator<field_t>(order));
	}

	/**
	 * \brief free memory allocated for the accelerators
	 */
	~field_type_accelerator_pool(void)
	{
		for (unsigned order = 0; order < MAX_ORDER; ++order)
			delete (*this)[order];
	}
};

#endif // FIELD_TYPE_ACCELERATOR_HPP_INCLUDED
