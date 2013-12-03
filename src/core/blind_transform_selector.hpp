#ifndef BLIND_TRANSFORM_SELECTOR_HPP_INCLUDED
#define BLIND_TRANSFORM_SELECTOR_HPP_INCLUDED

#include "singularity_types.hpp"
#include "domain.hpp"

namespace blind_transform
{
	struct duffy {};
	struct square {};
}

template <class singularity_type, class domain>
struct blind_transform_selector;

template <>
struct blind_transform_selector<singularity_type::log<1>, line_domain>
{
	typedef blind_transform::square type;
};


#endif // BLIND_TRANSFORM_SELECTOR_INCLUDED
