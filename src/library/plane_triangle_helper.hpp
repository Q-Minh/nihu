/** \file plane_triangle_helper.hpp
 * \brief helper functions to compute analytical integrals over plane triangles
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 */
 
#ifndef PLANE_TRIANGLE_HELPER_HPP_INCLUDED
#define PLANE_TRIANGLE_HELPER_HPP_INCLUDED

/**
 * \brief compute angles and radii in a plane triangle
 * \tparam scalar_t the scalar type
 * \param [in] elem the element
 * \param [out] r the radii from the collocation point
 * \param [out] theta the angles between the radii and the sides
 * \param [out] alpha the angles between the radii
 */
template <class scalar_t>
void planar_triangle_helper(
	tria_1_elem const &elem,
	scalar_t r[],
	scalar_t theta[],
	scalar_t alpha [])
{
	auto const &C_old = elem.get_coords();
	auto const &x0 = elem.get_center();
	unsigned const N = tria_1_elem::num_nodes;

	typename tria_1_elem::coords_t R, C;
	for (unsigned i = 0; i < N; ++i)
	{
		R.col(i) = C_old.col(i) - x0;
		r[i] = R.col(i).norm();
		R.col(i) /= r[i];
		C.col(i) = C_old.col(i) - C_old.col((i+1) % N);
		C.col(i) /= C.col(i).norm();
	}

	for (unsigned i = 0; i < N; ++i)
	{
		theta[i] = std::acos(R.col(i).dot(R.col((i+1) % N)));
		alpha[i] = std::acos(R.col(i).dot(C.col(i)));
	}
}

#endif // PLANE_TRIANGLE_HELPER_HPP_INCLUDED
