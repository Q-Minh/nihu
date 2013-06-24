#ifndef CUSTOM_GAUSSIAN_ELEMENT_HPP
#define CUSTOM_GAUSSIAN_ELEMENT_HPP

#include "shapeset.hpp"


class quad_1_gauss_shape_set;

template<>
struct shape_set_traits<quad_1_gauss_shape_set>
{
	typedef quad_domain domain_t;
	static unsigned const num_nodes = 4;
	static unsigned const polynomial_order = 1;
	static unsigned const jacobian_order = 1;
};


class quad_1_gauss_shape_set
	: public shape_set_base<quad_1_gauss_shape_set>
{
public:
	static shape_t eval_shape(xi_t const &_xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1];
		shape_t L;
		L <<
			(1.0 - sqrt(3.0)*xi) * (1.0 - sqrt(3.0)*eta) / 4.0,
			(1.0 + sqrt(3.0)*xi) * (1.0 - sqrt(3.0)*eta) / 4.0,
			(1.0 + sqrt(3.0)*xi) * (1.0 + sqrt(3.0)*eta) / 4.0,
			(1.0 - sqrt(3.0)*xi) * (1.0 + sqrt(3.0)*eta) / 4.0;
		return L;
	}

	static dshape_t eval_dshape(xi_t const & _xi)
	{
		scalar_t xi = _xi[0], eta = _xi[1];
		dshape_t dL;
		dL <<
			(-sqrt(3.0)) * (1.0 - sqrt(3.0)*eta) / 4.0, (1.0 - sqrt(3.0)*xi) * (-sqrt(3.0)) / 4.0,
			(+sqrt(3.0)) * (1.0 - sqrt(3.0)*eta) / 4.0, (1.0 + sqrt(3.0)*xi) * (-sqrt(3.0)) / 4.0,
			(+sqrt(3.0)) * (1.0 + sqrt(3.0)*eta) / 4.0, (1.0 + sqrt(3.0)*xi) * (+sqrt(3.0)) / 4.0,
			(-sqrt(3.0)) * (1.0 + sqrt(3.0)*eta) / 4.0, (1.0 - sqrt(3.0)*xi) * (+sqrt(3.0)) / 4.0;
		return dL;
	}

	static xi_t const *corner_begin(void)
	{
		return m_corners;
	}

protected:
	static const xi_t m_corners[num_nodes];
};


quad_1_gauss_shape_set::xi_t
	const quad_1_gauss_shape_set::m_corners[quad_1_gauss_shape_set::num_nodes] = {
		quad_1_gauss_shape_set::xi_t(-sqrt(3.0)/3.0,-sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(+sqrt(3.0)/3.0,-sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(+sqrt(3.0)/3.0,+sqrt(3.0)/3.0),
		quad_1_gauss_shape_set::xi_t(-sqrt(3.0)/3.0,+sqrt(3.0)/3.0)
};


#include "duffy_quadrature.hpp"

/** \brief Specialisation of ::duffy_traits for ::quad_1_gauss_shape_set */
template <>
struct duffy_traits<quad_1_gauss_shape_set>
{
	/** \brief indices of the Duffy corners for singular corners */
	static unsigned const duffy_corner_indices[4][3+1];
	/** \brief indices of the Duffy corners for internal singular point */
	static unsigned const duffy_face_indices[5+1];
};

unsigned const duffy_traits<quad_1_gauss_shape_set>::duffy_corner_indices[4][3+1] = {
	{2, /*|*/ 1, 2, 3},
	{2, /*|*/ 2, 3, 0},
	{2, /*|*/ 3, 0, 1},
	{2, /*|*/ 0, 1, 2}
};

unsigned const duffy_traits<quad_1_gauss_shape_set>::duffy_face_indices[5+1] = {
	4, /*|*/ 0, 1, 2, 3, 0
};

#include "element.hpp"

typedef general_surface_element<quad_1_gauss_shape_set, 1000> quad_1_gauss_elem;

#endif // CUSTOM_GAUSSIAN_ELEMENT_HPP

