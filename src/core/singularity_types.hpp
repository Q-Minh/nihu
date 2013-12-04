#ifndef SINGULARITY_TYPES_HPP_INCLUDED
#define SINGULARITY_TYPES_HPP_INCLUDED

namespace singularity_type
{
	struct regular {};

	template <unsigned order>
	struct log {};

	template <unsigned order>
	struct inverse {};
}

template <class SingularityType>
struct minimal_reference_dimension;

template <>
struct minimal_reference_dimension<singularity_type::log<1> >
{
	static unsigned const value = 1;
};


template <unsigned n>
struct minimal_reference_dimension<singularity_type::inverse<n> >
{
	static unsigned const value = n+1;
};


#endif // SINGULARITY_TYPES_HPP_INCLUDED
