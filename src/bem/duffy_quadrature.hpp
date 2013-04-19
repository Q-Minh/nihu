#ifndef DUFFY_QUADRATURE_HPP_INCLUDED
#define DUFFY_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"

template <class LSet>
struct duffy_traits;

template <>
struct duffy_traits<tria_1_shape_set>
{
	static unsigned const duffy_indices[3][3];
};

unsigned const duffy_traits<tria_1_shape_set>::duffy_indices[3][3] = {
	{1,  1, 2},
	{1,  2, 0},
	{1,  0, 1}
};


template <>
struct duffy_traits<quad_1_shape_set>
{
	static unsigned const duffy_indices[4][4];
};

unsigned const duffy_traits<quad_1_shape_set>::duffy_indices[4][4] = {
	{2, /*|*/ 1, 2, 3},
	{2, /*|*/ 2, 3, 0},
	{2, /*|*/ 3, 0, 1},
	{2, /*|*/ 0, 1, 2}
};

template <class LSet>
typename quadrature_type<gauss_family_tag, typename LSet::domain_t>::type
	duffy_quadrature_corner(unsigned degree, unsigned singular_corner)
{
	typedef typename LSet::domain_t domain_t;
	typedef typename quadrature_type<gauss_family_tag, domain_t>::type ret_quad_type;
	typedef typename quadrature_type<gauss_family_tag, quad_domain>::type source_quad_type;
	typedef Eigen::Matrix<typename domain_t::scalar_t, 4, domain_t::dimension> coords_t;

	unsigned const *array = duffy_traits<LSet>::duffy_indices[singular_corner];
	unsigned num_duffies = *array;
	array++;

	ret_quad_type result;
	source_quad_type source(degree);

	coords_t coords;
	coords.row(0) = coords.row(1) = LSet::corner_at(singular_corner);

	for (size_t d = 0; d < num_duffies; ++d)
	{
		coords.row(2) = LSet::corner_at(array[d]);
		coords.row(3) = LSet::corner_at(array[d+1]);

		result += source.transform<quad_1_shape_set>(coords);
	}

	return result;
}

#endif // DUFFY_QUADRATURE_HPP_INCLUDED
