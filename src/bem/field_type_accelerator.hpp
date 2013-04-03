#ifndef FIELD_TYPE_ACCELERATOR_HPP_INCLUDED
#define FIELD_TYPE_ACCELERATOR_HPP_INCLUDED

#include "field.hpp"
#include "quadrature.hpp"

#include <Eigen/StdVector>
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >

template <class Field>
class field_type_accelerator
{
public:
	typedef Field field_t;
	typedef typename field_t::elem_t elem_t;
	typedef typename elem_t::domain_t domain_t;
	typedef gauss_quadrature<domain_t> quadrature_t;
	typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;

	typedef typename field_t::nset_t nset_t;
	typedef typename nset_t::shape_t shape_t;
	typedef EIGENSTDVECTOR(shape_t) nset_vector_t;

	class field_accelerator_element
	{
	public:
		field_accelerator_element(quadrature_elem_t const *qe, shape_t const *sh) : m_quadrature_element(qe), m_shape(sh)
		{
		}

		quadrature_elem_t const &get_quadrature_elem(void) const
		{
			return *m_quadrature_element;
		}

		shape_t const &get_shape(void) const
		{
			return *m_shape;
		}

	private:
		quadrature_elem_t const *m_quadrature_element;
		shape_t const *m_shape;
	};

	class iterator
	{
	public:
		iterator (field_type_accelerator const &parent, unsigned i)
			: m_parent(parent), m_idx(i), m_elem(&m_parent.m_quadrature[m_idx], &m_parent.m_nset_vector[m_idx])
		{
		}

		iterator &operator++()
		{
			++m_idx;
			m_elem = field_accelerator_element(&m_parent.m_quadrature[m_idx], &m_parent.m_nset_vector[m_idx]);
			return *this;
		}

		bool operator==(iterator const &other)
		{
			return &m_parent == &other.m_parent && m_idx == other.m_idx;
		}

		bool operator!=(iterator const &other)
		{
			return ! (*this == other);
		}

		field_accelerator_element const &operator*() const
		{
			return m_elem;
		}

		field_accelerator_element const *operator->() const
		{
			return &m_elem;
		}

	private:
		field_type_accelerator const &m_parent;
		unsigned m_idx;
		field_accelerator_element m_elem;
	};

	iterator begin(void) const
	{
		return iterator(*this, 0);
	}

	iterator end(void) const
	{
		return iterator(*this, m_quadrature.size());
	}

	field_type_accelerator(unsigned order) : m_quadrature(order)
	{
		m_nset_vector.reserve(m_quadrature.size());
		for (auto quad_it = m_quadrature.begin(); quad_it != m_quadrature.end(); ++quad_it)
			m_nset_vector.push_back(nset_t::eval_shape(quad_it->get_xi()));
	}

protected:
	quadrature_t m_quadrature;
	nset_vector_t m_nset_vector;
};

#endif // FIELD_TYPE_ACCELERATOR_HPP_INCLUDED
