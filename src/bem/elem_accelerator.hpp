#ifndef ELEM_ACCELERATOR_HPP
#define ELEM_ACCELERATOR_HPP

#include <algorithm>
#include <iterator>

#include <Eigen/StdVector>
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >

#include "quadrature.hpp"
#include "descriptor.hpp"

template <class ElemDescriptor>
class elem_accelerator : public EIGENSTDVECTOR(ElemDescriptor)
{
public:
	typedef ElemDescriptor elem_descriptor_t;
	
	template <class elem_t>
	elem_accelerator(elem_t const &elem, quadrature<typename elem_t::domain_t> const &q)
	{
		typedef quadrature<typename elem_t::domain_t> quadrature_t;
		typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;
		
		this->reserve(q.size());
		std::transform(q.begin(), q.end(),
			std::back_inserter(*this),
			[&elem] (quadrature_elem_t const &qe) { return elem_descriptor_t(elem, qe); }
		);
	}
};

template <class ElemType, class KernelInput>
class elem_pool : public std::vector<elem_accelerator<KernelInput> >
{
public:
	typedef ElemType elem_t;
	typedef KernelInput kernel_input_t;

	typedef gauss_quadrature<typename elem_t::domain_t> quadrature_t;

	template <class quadIter>
	elem_pool(elem_t const &elem, quadIter begin, quadIter end)
	{
		std::transform(begin, end,
			std::back_inserter(*this),
			[&elem] (quadrature_t const &q) { return elem_accelerator<kernel_input_t>(elem, q); }
		);
	}
};

#endif

