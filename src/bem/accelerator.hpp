#ifndef ACCELERATOR_HPP_INCLUDED
#define ACCELERATOR_HPP_INCLUDED

#include <vector>

#include "elem_accelerator.hpp"
#include "nset_accelerator.hpp"

template <class ElemType, class  KernelInput>
class accelerator
{
public:
	typedef ElemType elem_t;
	typedef KernelInput kernel_input_t;

	typedef typename elem_t::domain_t domain_t;
	typedef gauss_quadrature<domain_t> quadrature_t;
	typedef typename elem_t::lset_t nset_t;

	typedef elem_pool<elem_t, kernel_input_t> elem_pool_t;
	typedef nset_pool<nset_t> nset_pool_t;

	accelerator()
	{
		std::vector<quadrature_t> q;
		q.push_back(quadrature_t(1));
		q.push_back(quadrature_t(2));
		q.push_back(quadrature_t(3));
		n_pool = nset_pool_t(q.begin(), q.end());
	}

protected:
	std::vector<elem_pool_t> elem_pools;
	nset_pool_t n_pool;
};

#endif

