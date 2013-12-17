#ifndef GUIGGIANI_1992_HPP_INCLUDED
#define GUIGGIANI_1992_HPP_INCLUDED

template <class Kernel, class Field>
class guiggiani_hypersingular_collocation
{
	typedef Field field_t;

	typedef typename field_t::elem_t elem_t;
	typedef typename field_t::nset_t nset_t;

	typedef typename nset_t::shape_t n_shape_t;

	typedef typename elem_t::xi_t xi_t;

	typename elem_t::x_t A_vector(elem_t const &elem, xi_t const &xi0, double theta)
	{
		auto dx = elem.get_dx(xi0);
		return dx.col(0) * cos(theta) + dx.col(1) * sin(theta);
	}

	n_shape_t N0(xi_t const &xi0)
	{
		return nset_t::eval_shape(xi0);
	}

	n_shape_t N1(xi_t const &xi0, double theta)
	{
		auto dN = nset_t::eval_dshape(xi0);
		return dN.col(0)*cos(theta) + dN.col(1)*sin(theta);
	}
};

#endif // GUIGGIANI_1992_HPP_INCLUDED
