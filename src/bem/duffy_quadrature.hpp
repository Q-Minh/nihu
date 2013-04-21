/**
 * \file duffy_quadrature.hpp
 * \brief Duffy-type singular quadrature transformations
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef DUFFY_QUADRATURE_HPP_INCLUDED
#define DUFFY_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"

// forward declaration
template <class LSet>
struct duffy_traits;


/**
 * \brief Transform regular quadratures into weakly singular Duffy-type quadratures
 * \tparam QuadFamily the regular quadrature family
 * \tparam LSet the element geometrical shape set representation
 * \todo the on_face and on_corner functions are almost identical. Unify them!
 */
template <class QuadFamily, class LSet>
class duffy_quadrature
{
	// CRTP check
	static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
		"The LSet parameter must be derived from shape_set_base<LSet>");
public:
	/** \brief template argument as nested type */
	typedef QuadFamily quadrature_family_tag;
	/** \brief template argument as nested type */
	typedef LSet lset_t;

	/** \brief element's domain type */
	typedef typename lset_t::domain_t domain_t;
	/** \brief base domain variable type */
	typedef typename domain_t::xi_t xi_t;
	/** \brief domain's scalar type */
	typedef typename domain_t::scalar_t scalar_t;

	/** \brief the return quadrature type */
	typedef typename quadrature_type<quadrature_family_tag, domain_t>::type ret_quad_type;
	/** \brief the source quadrature type (regular quad quadrature) */
	typedef typename quadrature_type<quadrature_family_tag, quad_domain>::type source_quad_type;

	/** \brief the local coordinate matrix type of the transformation */
	typedef Eigen::Matrix<scalar_t, 4, domain_t::dimension> coords_t;

	/**
	 * \brief create a duffy quadrature that is singular on one corner of the selected element
	 * \param [in] degree the polynomial degree of the original regular quadrature
	 * \param [in] singular_corner index of the singular corner
	 * \return a Duffy type singular quadrature
	 */
	static ret_quad_type on_corner(unsigned degree, unsigned singular_corner)
	{
		unsigned const *array = duffy_traits<LSet>::duffy_corner_indices[singular_corner];
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

			result += source.template transform<quad_1_shape_set>(coords);
		}

		return result;
	}

	/**
	 * \brief create a duffy quadrature that is singular on the face of the selected element
	 * \param [in] degree the polynomial degree of the original regular quadrature
	 * \param [in] singular_point coordinates of the singular point
	 * \return a Duffy type singular quadrature
	 */
	static ret_quad_type on_face(unsigned degree, xi_t const &singular_point)
	{
		unsigned const *array = duffy_traits<LSet>::duffy_face_indices;
		unsigned num_duffies = *array;
		array++;

		ret_quad_type result;
		source_quad_type source(degree);

		coords_t coords;
		coords.row(0) = coords.row(1) = singular_point;

		for (size_t d = 0; d < num_duffies; ++d)
		{
			coords.row(2) = LSet::corner_at(array[d]);
			coords.row(3) = LSet::corner_at(array[d+1]);

			result += source.template transform<quad_1_shape_set>(coords);
		}

		return result;
	}
}; // end of class duffy_quadrature



template <>
struct duffy_traits<tria_1_shape_set>
{
	static unsigned const duffy_corner_indices[3][2+1];
	static unsigned const duffy_face_indices[4+1];
};

unsigned const duffy_traits<tria_1_shape_set>::duffy_corner_indices[3][2+1] = {
	{1, /*|*/ 1, 2},
	{1, /*|*/ 2, 0},
	{1, /*|*/ 0, 1}
};

unsigned const duffy_traits<tria_1_shape_set>::duffy_face_indices[4+1] = {
	3, /*|*/ 0, 1, 2, 0
};


template <>
struct duffy_traits<quad_1_shape_set>
{
	static unsigned const duffy_corner_indices[4][3+1];
	static unsigned const duffy_face_indices[5+1];
};

unsigned const duffy_traits<quad_1_shape_set>::duffy_corner_indices[4][3+1] = {
	{2, /*|*/ 1, 2, 3},
	{2, /*|*/ 2, 3, 0},
	{2, /*|*/ 3, 0, 1},
	{2, /*|*/ 0, 1, 2}
};

unsigned const duffy_traits<quad_1_shape_set>::duffy_face_indices[5+1] = {
	4, /*|*/ 0, 1, 2, 3, 0
};


#endif // DUFFY_QUADRATURE_HPP_INCLUDED
