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
	typedef NSet nset_t;

	typedef typename nset_t::domain_t domain_t;
	typedef gauss_quadrature<domain_t> quadrature_t;
	typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;

	nset_accelerator(quadrature_t const &q)
	{
		this->reserve(q.size());
		std::transform(
			q.begin(), q.end(),
			std::back_inserter(*this),
			[] (quadrature_elem_t const &qe) { return nset_t::eval_L(qe.get_xi()); }
		);
	}
};

template <class NSet>
class nset_pool : public std::vector<nset_accelerator<NSet> >
{
public:
	typedef NSet nset_t;
	typedef gauss_quadrature<typename nset_t::domain_t> quadrature_t;

	nset_pool()
	{
	}

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

