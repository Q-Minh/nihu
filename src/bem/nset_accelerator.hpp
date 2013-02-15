#ifndef NSET_ACCELERATOR_HPP
#define NSET_ACCELERATOR_HPP

#include <algorithm>
#include <iterator>

#include <Eigen/StdVector>
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >

#include "quadrature.hpp"
#include "shapeset.hpp"

template <class NSet>
class nset_accelerator : public EIGENSTDVECTOR(typename NSet::L_t)
{
public:
	typedef NSet nset_t;	/**< template parameter as nested type */

	typedef typename nset_t::domain_t domain_t; /**< \brief the n-set's domain type */
	typedef gauss_quadrature<domain_t> quadrature_t;	/**< \brief quadrature type of the domain */
	typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;	/**< \brief quadrature element type */

	/**
	 * \brief construct from a quadrature
	 * \param q a quadrature object
	 */
	nset_accelerator(quadrature_t const &q)
	{
		this->reserve(q.size());
		std::transform(
			q.begin(), q.end(),
			std::back_inserter(*this),
			[] (quadrature_elem_t const &qe) { return nset_t::eval_shape(qe.get_xi()); }
		);
	}
};

template <class NSet>
class nset_pool : public std::vector<nset_accelerator<NSet> >
{
public:
	typedef NSet nset_t;	/**< \brief template parameter as nested type */
	typedef typename nset_t::domain_t domain_t; /**< \brief the n-set's domain type */
	typedef gauss_quadrature<domain_t> quadrature_t;	/**< \brief quadrature type of the domain */

	nset_pool()
	{
	}

	/**
	 * \brief construct from a set of quadratures
	 * \param begin begin iterator of the quadratures
	 * \param end end iterator of quadratures
	 */
	template <class quadIter>
	nset_pool(quadIter begin, quadIter end)
	{
		std::transform(begin, end,
			std::back_inserter(*this),
			[] (quadrature_t const &q) { return nset_accelerator<nset_t>(q); }
		);
	}
};


#endif

