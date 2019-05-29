/** \file nd_cheb.hpp
 * \brief n-dimensional Chebyshev polynomials
 */
#ifndef ND_CHEB_H_INCLUDED
#define ND_CHEB_H_INCLUDED

#include <boost/math/special_functions/chebyshev.hpp>

#include "bounding_box.hpp"

#include "util/math_functions.hpp"
#include "util/misc.hpp"


#include <Eigen/Dense>

namespace NiHu
{
namespace fmm
{

/** \brief linear transform of points from a bounding box to an other
 * \tparam Dim dimension of space
 * \tparam T the scalar type
 * \param [in] x Dim x N matrix of points to transform
 * \param [in] bbTo the destination bounding box
 * \param [in] bbFrom the source bounding box
 */
template <size_t Dim, class T = double>
Eigen::Matrix<T, Dim, Eigen::Dynamic> lintrans(
	Eigen::Matrix<T, Dim, Eigen::Dynamic> const &x,
	bounding_box<Dim> const &bbTo,
	bounding_box<Dim> const &bbFrom)
{
	Eigen::Matrix<T, Dim, Eigen::Dynamic> y = x;
	for (size_t d = 0; d < Dim; ++d)
	{
		y.row(d).array() -= bbFrom.get_center()(d);
		y.row(d).array() *= bbTo.get_diameter() / bbFrom.get_diameter();
		y.row(d).array() += bbTo.get_center()(d);
	}
	return y;
}

/**
 * \brief return Chebyshev nodes in a multidimensional bounding box
 * \tparam Dim space dimension
 * \tparam T the scalar type
 * \param [in] n required number of Chebyshev nodes
 * \param [in] bb the bounding box of Chebyshev nodes
 */
template <size_t Dim, class T = double>
Eigen::Matrix<T, Dim, Eigen::Dynamic> chebnodes(size_t n,
	bounding_box<Dim> const &bb = bounding_box<Dim>())
{
	// get Chebyshev roots in [-1 +1]
	Eigen::Matrix<T, 1, Eigen::Dynamic> x0(1, n);
	for (size_t i = 0; i < n; ++i)
		x0(i) = chebroot(n, i);

	// nd Chebyshev nodes in [-1 +1]^D
	Eigen::Matrix<T, Dim, Eigen::Dynamic> x(Dim, Ipow(n, Dim));
	for (size_t i = 0; i < Ipow(n, Dim); ++i)
	{
		size_t j = i;
		for (size_t d = 0; d < Dim; ++d)
		{
			x(d, i) = x0(j % n);
			j /= n;
		}
	}

	// transform to bb
	x = lintrans<Dim, T>(x, bb, bounding_box<Dim>());
	return x;
}

/**
 * \brief 1D Chebyshev anterpolation matrix
 * \tparam T the scalar type
 * \param [in] N Chebyshev order
 * \param [in] y source nodes of anterpolation
 * \param [in] x receiver nodes of anterpolation
 */
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
chebanterp0(size_t N, Eigen::Matrix<T, Eigen::Dynamic, 1> const &x,
	Eigen::Matrix<T, Eigen::Dynamic, 1> const &y)
{
	using boost::math::chebyshev_t;
	size_t n = x.rows(), m = y.rows();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> S(n, m);
	Eigen::Matrix<T, Eigen::Dynamic, 1> Sx(n, 1), Sy(m, 1);
	S.setConstant(1.);
	for (unsigned l = 1; l < N; ++l)
	{
		for (size_t i = 0; i < n; ++i)
			Sx(i) = chebyshev_t(l, x(i));
		for (size_t j = 0; j < m; ++j)
			Sy(j) = chebyshev_t(l, y(j));
		S += 2. * (Sx * Sy.transpose());
	}
	return S / N;
}

/**
 * \brief y-derivative 1D Chebyshev anterpolation matrix
 * \tparam T the scalar type
 * \param [in] N Chebyshev order
 * \param [in] y source nodes of anterpolation
 * \param [in] x receiver nodes of anterpolation
 */
template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
chebanterp0_dy(size_t N, Eigen::Matrix<T, Eigen::Dynamic, 1> const &x,
	Eigen::Matrix<T, Eigen::Dynamic, 1> const &y)
{
	using boost::math::chebyshev_t;
	using boost::math::chebyshev_t_prime;
	size_t n = x.rows(), m = y.rows();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> S(n, m);
	Eigen::Matrix<T, Eigen::Dynamic, 1> Sx(n, 1), Sy(m, 1);
	S.setZero();
	for (unsigned l = 1; l < N; ++l)
	{
		for (size_t i = 0; i < n; ++i)
			Sx(i) = chebyshev_t(l, x(i));
		for (size_t j = 0; j < m; ++j)
			Sy(j) = chebyshev_t_prime(l, y(j));
		S += 2. * (Sx * Sy.transpose());
	}
	return S / N;
}


/**
 * \brief ND Chebyshev anterpolation matrix from points to Chebyshev nodes
 * \tparam T the scalar type
 * \tparam Dim bounding box dimension
 * \param [in] N Chebyshev order
 * \param [in] bb the bounding box of Chebyshev nodes
 * \param [in] y0 the source nodes of anterpolation
 */
template <class T, size_t Dim>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
chebanterp(size_t N, bounding_box<Dim> const &bb,
	Eigen::Matrix<T, Dim, Eigen::Dynamic> const &y0)
{
	auto x = chebnodes<Dim, T >(N);
	auto y = lintrans<Dim, T>(y0, bounding_box<Dim>(), bb);
	size_t n = x.cols(), m = y.cols();
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> S(n, m);
	S.setConstant(1.);
	for (size_t d = 0; d < Dim; ++d)
		S.array() *= chebanterp0<T>(N, x.row(d).transpose(), y.row(d).transpose()).array();
	return S;
}

/**
 * \brief normal derivative of Chebyshev anterpolation matrix
 * \tparam T the scalar type
 * \tparam Dim bounding box dimension
 * \param [in] N Chebyshev order
 * \param [in] bb the bounding box of Chebyshev nodes
 * \param [in] y0 the source nodes of anterpolation
 * \param [in] n0 the source normal vectors
 */
template <class T, size_t Dim>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
chebanterp_dny(size_t N,
	bounding_box<Dim> const &bb,
	Eigen::Matrix<T, Dim, Eigen::Dynamic> const &y0,
	Eigen::Matrix<T, Dim, Eigen::Dynamic> const &ny0
)
{
	auto x = chebnodes<Dim, T >(N);
	auto y = lintrans<Dim, T>(y0, bounding_box<Dim>(), bb);
	size_t n = x.cols(), m = y.cols();
	double scale = bb.get_diameter() / bounding_box<Dim>().get_diameter();

	Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> dS(n, m);
	dS.setZero();

	for (size_t q = 0; q < Dim; ++q)
	{
		Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> S(n, m);
		S.setConstant(1.);
		for (size_t d = 0; d < Dim; ++d)
		{
			if (d == q)
				S *= chebanterp0_dy<T>(N, x.row(d).transpose(), y.row(d).transpose()).array();
			else
				S *= chebanterp0<T>(N, x.row(d).transpose(), y.row(d).transpose()).array();
		}
		dS += S.rowwise() * ny0.row(q).array();
	}
	return dS / scale;
}

} // end of namespace fmm
} // namespace NiHu

#endif //ND_CHEB_H_INCLUDED
