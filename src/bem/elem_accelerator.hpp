#ifndef ELEM_ACCELERATOR_HPP
#define ELEM_ACCELERATOR_HPP

#include <algorithm>

#include "quadrature.hpp"
#include "element.hpp"
#include "descriptor.hpp"

template <ElemDescriptor>
class elem_accelerator : public EIGENSTDVECTOR(ElemDescriptor)
{
public:
	typedef ElemDescriptor elem_descriptor_t
	
	template <class elem_t>
	elem_accelerator(elem_t const &elem, quadrature<elem_t::domain_t> const &quad)
	{
		typedef quadrature<elem_t::domain_t> quadrature_t;
		typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;
		
		this->resize(quad.size());
		std::for_each(
			quad.begin(),
			quad.end(),
			[this, &elem] (quadrature_elem_t const &qe) {
				this->push_back(elem_descriptor_t(elem, qe.get_xi()));
			}
		);
	}

protected:
};

#endif

