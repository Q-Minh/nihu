/**
 * \file quadrature.hpp
 * \author Peter fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 * \brief Decalaration of class gauss_quad and its specialisations
 */
#ifndef QUADRATURE_HPP_INCLUDED
#define QUADRATURE_HPP_INCLUDED

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include "domain.hpp"

template <class Domain, unsigned N>
class gauss_quad;

template <unsigned N>
class gauss_quad<line_domain, N>
{
public:
	static const unsigned size = N;
	typedef Eigen::Matrix<double, size, 1> xivec_t;
	typedef Eigen::Matrix<double, size, 1> weightvec_t;

	static void init(void)
	{
		typedef Eigen::Matrix<double, size, size> Mat_t;
		Mat_t A = Mat_t::Zero();
		for (unsigned i = 1; i < size; ++i)
			A(i, i-1) = A(i-1, i) = i / sqrt(4.0*(i*i)-1.0);
		Eigen::SelfAdjointEigenSolver<Mat_t> S(A);
		xi = S.eigenvalues();
		w = 2.0 * S.eigenvectors().row(0).cwiseAbs2().transpose();
	}
	
	static xivec_t const &get_xi(void)
	{
		return xi;
	}
	
	static weightvec_t const &get_weight(void)
	{
		return w;
	}

protected:
	static xivec_t xi;
	static weightvec_t w;
};

template <unsigned N>
typename gauss_quad<line_domain, N>::xivec_t gauss_quad<line_domain, N>::xi;
template <unsigned N>
typename gauss_quad<line_domain, N>::weightvec_t gauss_quad<line_domain, N>::w;


template <unsigned N>
class gauss_quad<quad_domain, N>
{
public:
	static const unsigned size = N*N;
	typedef Eigen::Matrix<double, size, 2> xivec_t;
	typedef Eigen::Matrix<double, size, 1> weightvec_t;

	static void init(void)
	{
		gauss_quad<line_domain, N>::init();
		
		typename gauss_quad<line_domain, N>::xivec_t _xi = gauss_quad<line_domain, N>::get_xi();
		typename gauss_quad<line_domain, N>::weightvec_t _w = gauss_quad<line_domain, N>::get_weight();

		typedef typename xivec_t::Index index_t;
		for (index_t i = 0, k = 0; i < N; ++i)
		{
			for (index_t j = 0; j < N; ++j, ++k)
			{
				xi(k,0) = _xi(i);
				xi(k,1) = _xi(j);
				w(k) = _w(i)*_w(j);
			}
		}
	}
	
	static xivec_t const &get_xi(void)
	{
		return xi;
	}
	
	static weightvec_t const &get_weight(void)
	{
		return w;
	}

protected:
	static xivec_t xi;
	static weightvec_t w;
};

template <unsigned N>
typename gauss_quad<quad_domain, N>::xivec_t gauss_quad<quad_domain, N>::xi;
template <unsigned N>
typename gauss_quad<quad_domain, N>::weightvec_t gauss_quad<quad_domain, N>::w;

#endif

