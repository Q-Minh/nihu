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
/**
 * \brief mapping from physical to intrinsic coordinates
 * \tparam Elem the element type
 */
template <class Elem>
class inverse_mapping;

/** \brief inverse mapping for surface elements
 * \tparam LSet the geometrical shape set
 * \tparam Scalar the element's scalar type
 */
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
	/** \brief constructor
	 * \param [in] elem the element
	 */
	inverse_mapping(elem_t const &elem)
		: m_elem(elem)
	{
	}

	/** \brief compute inverse mapping
	 * \param [in] x0 physical coordinates
	 * \param [in] tol tolerance (absolute error in physical coordinates)
	 * \param [in] max_iter maximal number of iterations
	*/
	bool eval(x_t const &x0, double tol, unsigned max_iter)
	{
		// initial guess
		m_res.setZero();
		
		// Newton Raphson iterations
		for (m_iter = 0; m_iter < max_iter; ++m_iter)
		{
			// intrinsic coordinates within the element
			auto xi = m_res.topRows(xi_dim);
			double zeta = m_res(xi_dim);

			// compute physical coordinates of guess
			auto y = m_elem.get_x(xi);
			x_t n = m_elem.get_normal(xi);
			auto x_guess = y + n * zeta;

			// check if converged
			auto f = x_guess - x0;
			m_error = f.norm();
			if (m_error < tol)
				return true;

			// compute derivatives
			auto dy = m_elem.get_dx(xi);
			auto ddy = m_elem.get_ddx(xi);

			// Newton Raphson derivative matrix
			Eigen::Matrix<double, N, N> df;

			if (N == 3)
			{
				auto const &dyxi = dy.col(shape_derivative_index::dXI);
				auto const &dyeta = dy.col(shape_derivative_index::dETA);

				auto const &ddyxixi = ddy.col(shape_derivative_index::dXIXI);
				auto const &ddyxieta = ddy.col(shape_derivative_index::dXIETA);
				auto const &ddyetaeta = ddy.col(shape_derivative_index::dETAETA);

				x_t dnxi = ddyxixi.cross(dyeta) + dyxi.cross(ddyxieta);
				x_t dneta = ddyxieta.cross(dyeta) + dyxi.cross(ddyetaeta);

				df.col(0) = dyxi + dnxi * zeta;
				df.col(1) = dyeta + dneta * zeta;
				df.col(2) = n;
			}
			else
				throw std::logic_error("Inverse mapping is unimplemented fr this element type");

			m_res = m_res - df.colPivHouseholderQr().solve(f);
		}
		
		return false;
	}

	/** \brief get the result of inverse mapping */
	result_t const &get_result() const
	{
		return m_res;
	}

	/** \brief get the number of iterations */
	unsigned get_iter() const
	{
		return m_iter;
	}

	/** \brief get the error (absolute in physical coordinates) */
	double get_error() const
	{
		return m_error;
	}

private:
	elem_t m_elem;
	result_t m_res;
	unsigned m_iter;
	double m_error;
};

} // end of namespace NiHu
 
#endif // INVERSE_MAPPING_HPP_INCLUDED
