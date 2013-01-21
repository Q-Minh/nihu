/**
 * \file quadrature.hpp
 * \author Peter fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 * \brief Decalaration of class quadrature, class gauss_quad and its specialisations
 */
#ifndef QUADRATURE_HPP_INCLUDED
#define QUADRATURE_HPP_INCLUDED

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "domain.hpp"

/**
 * \brief store quadrature base points and weights
 * \tparam Domain the base domain where the quadrature points are defined
 * \tparam Size number of quadrature points
 */
template <class Domain, unsigned Size>
class quadrature
{
public:
    /** \brief template parameter as nested type */
	typedef Domain domain_t;
    /** \brief template parameter as nested constant */
	static unsigned const size = Size;

    /** \brief scalar type of the domain */
	typedef typename domain_t::scalar_t scalar_t;

    /** \brief type of the base points matrix */
	typedef Eigen::Matrix<scalar_t, size, domain_t::dimension> xivec_t;
    /** \brief type of the weights vector */
	typedef Eigen::Matrix<scalar_t, size, 1> weightvec_t;

    /**
      * \brief return base points
      * \return reference to statically stored base points
      */
	static xivec_t const &get_xi(void)
	{
		return xi;
	}

    /**
      * \brief return weights
      * \return reference to statically stored weights
      */
	static weightvec_t const &get_weight(void)
	{
		return w;
	}

protected:
	static xivec_t xi;
	static weightvec_t w;
};

template <class Domain, unsigned Size>
typename quadrature<Domain, Size>::xivec_t quadrature<Domain, Size>::xi;
template <class Domain, unsigned Size>
typename quadrature<Domain, Size>::weightvec_t quadrature<Domain, Size>::w;


template <class Domain, unsigned N>
class gauss_quad;

template <unsigned N>
class gauss_quad<line_domain, N> : public quadrature<line_domain, N>
{
public:
	typedef quadrature<line_domain, N> base;

	static void init(void)
	{
		typedef Eigen::Matrix<typename base::scalar_t, base::size, base::size> Mat_t;
		Mat_t A = Mat_t::Zero();
		for (unsigned i = 1; i < base::size; ++i)
			A(i, i-1) = A(i-1, i) = i / sqrt(4.0*(i*i)-1.0);
		Eigen::SelfAdjointEigenSolver<Mat_t> S(A);
		base::xi = S.eigenvalues();
		base::w = 2.0 * S.eigenvectors().row(0).cwiseAbs2().transpose();
	}
};

template <unsigned N>
class gauss_quad<quad_domain, N> : public quadrature<quad_domain, N*N>
{
public:
	typedef quadrature<quad_domain, N*N> base;

	static void init(void)
	{
		gauss_quad<line_domain, N>::init();

		typename gauss_quad<line_domain, N>::xivec_t _xi = gauss_quad<line_domain, N>::get_xi();
		typename gauss_quad<line_domain, N>::weightvec_t _w = gauss_quad<line_domain, N>::get_weight();

		typedef typename base::xivec_t::Index index_t;
		for (index_t i = 0, k = 0; i < N; ++i)
		{
			for (index_t j = 0; j < N; ++j, ++k)
			{
				base::xi(k,0) = _xi(i);
				base::xi(k,1) = _xi(j);
				base::w(k) = _w(i)*_w(j);
			}
		}
	}
};

template <>
class gauss_quad<tria_domain, 1> : public quadrature<tria_domain, 1>
{
public:
	typedef quadrature<tria_domain, 1> base;

	static void init(void)
	{
		base::xi <<
			1.0/3.0, 1.0/3.0;
		base::w <<
			1.0/2.0;
	}
};


template <>
class gauss_quad<tria_domain, 2> : public quadrature<tria_domain, 3>
{
public:
	typedef quadrature<tria_domain, 3> base;

	static void init(void)
	{
		base::xi <<
			1.0/6.0, 4.0/6.0,
			1.0/6.0, 1.0/6.0,
			4.0/6.0, 1.0/6.0;
		base::w <<
			1.0/6.0,
			1.0/6.0,
			1.0/6.0;
	}
};


template <>
class gauss_quad<tria_domain, 3> : public quadrature<tria_domain, 4>
{
public:
	typedef quadrature<tria_domain, 4> base;

	static void init(void)
	{
		base::xi <<
			1.0/3.0, 1.0/3.0,
			1.0/5.0, 3.0/5.0,
			1.0/5.0, 1.0/5.0,
			3.0/5.0, 1.0/5.0;
		base::w <<
			-0.281250000000000,
			 0.260416666666667,
			 0.260416666666667,
			 0.260416666666667;
	}
};


template <>
class gauss_quad<tria_domain, 4> : public quadrature<tria_domain, 6>
{
public:
	typedef quadrature<tria_domain, 6> base;

	static void init(void)
	{
		base::xi <<
			0.445948490915965, 0.108103018168070,
			0.445948490915965, 0.445948490915965,
			0.108103018168070, 0.445948490915965,
			0.091576213509771, 0.816847572980459,
			0.091576213509771, 0.091576213509771,
			0.816847572980459, 0.091576213509771;
		base::w <<
			0.111690794839006,
			0.111690794839006,
			0.111690794839006,
			0.054975871827661,
			0.054975871827661,
			0.054975871827661;
	}
};


template <>
class gauss_quad<tria_domain, 5> : public quadrature<tria_domain, 7>
{
public:
	typedef quadrature<tria_domain, 7> base;

	static void init(void)
	{
		base::xi <<
			0.333333333333333, 0.333333333333333,
			0.470142064105115, 0.059715871789770,
			0.470142064105115, 0.470142064105115,
			0.059715871789770, 0.470142064105115,
			0.101286507323456, 0.797426985353087,
			0.101286507323456, 0.101286507323456,
			0.797426985353087, 0.101286507323456;
		base::w <<
			0.112500000000000,
			0.066197076394253,
			0.066197076394253,
			0.066197076394253,
			0.062969590272414,
			0.062969590272414,
			0.062969590272414;
	}
};


#endif

