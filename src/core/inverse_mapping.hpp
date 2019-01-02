/**
 * \file inverse_mapping.hpp
 * \ingroup funcspace
 */

#ifndef INVERSE_MAPPING_HPP_INCLUDED
#define INVERSE_MAPPING_HPP_INCLUDED

#include <Eigen/Dense>

#include "element.hpp"
#include "shapeset.hpp"

namespace NiHu
{
template <class Elem>
class inverse_mapping;

template <class LSet, class Scalar>
class inverse_mapping<surface_element<LSet, Scalar> >
{
	typedef surface_element<LSet, Scalar> elem_t;
	typedef typename elem_t::xi_t xi_t;
	typedef typename elem_t::x_t x_t;
	static unsigned const xi_dim = elem_t::domain_t::dimension;
	static unsigned const N = xi_dim + 1;
	typedef Eigen::Matrix<Scalar, N, 1> result_t;
	
public:
	static result_t eval(elem_t const &elem, x_t const &x0, unsigned max_iter, double tol)
	{
		result_t res = result_t::Zero();
		
		for (unsigned i = 0; i < max_iter; ++i)
		{
			xi_t xi = res.topRows(xi_dim);
			double zeta = res(xi_dim);

			auto y = elem.get_x(xi);
			x_t n = elem.get_normal(xi);

			auto f = y + n * zeta - x0;
			if (f.norm() < tol)
				break;

			auto dy = elem.get_dx(xi);
			auto ddy = elem.get_ddx(xi);

			if (N == 3)
			{
				auto const &dyxi = dy.col(shape_derivative_index::dXI);
				auto const &dyeta = dy.col(shape_derivative_index::dETA);

				auto const &ddyxixi = ddy.col(shape_derivative_index::dXIXI);
				auto const &ddyxieta = ddy.col(shape_derivative_index::dXIETA);
				auto const &ddyetaeta = ddy.col(shape_derivative_index::dETAETA);

				x_t dnxi = ddyxixi.cross(dyeta) + dyxi.cross(ddyxieta);
				x_t dneta = ddyxieta.cross(dyeta) + dyxi.cross(ddyetaeta);

				Eigen::Matrix<double, 3, 3> df;
				df.col(0) = dyxi + dnxi * zeta;
				df.col(1) = dyeta + dneta * zeta;
				df.col(2) = n;

				res = res - df.colPivHouseholderQr().solve(f);
			}
		}
		
		return res;
	}
};

} // end of namespace NiHu

 
#endif // INVERSE_MAPPING_HPP_INCLUDED