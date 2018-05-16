// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file gaussian_quadrature.hpp
 * \ingroup quadrature
 * \brief implementation of Gaussian quadratures
 */

#ifndef GAUSSIAN_QUADRATURE_HPP_INCLUDED
#define GAUSSIAN_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"
#include "../library/lib_domain.hpp"

#include <stdexcept>

namespace NiHu
{

/**
 * \brief return 1D N-point Guassian quadrature
 * \tparam scalar_t the scalar type
 * \param [in] N number of Gaussian points
 * \return matrix containing the Gaussian locations and weights
 * \details The Gaussian locations are roots of the Legendre polynomials \f$\varphi_n(x)\f$.
 * The polynomials obey the recurrence relation
 *
 * \f$\varphi_{0}(x) = 1\f$
 *
 * \f$\varphi_{1}(x) = x\f$
 *
 * \f$\varphi_{n+1}(x) = x\varphi_{n}(x) - b_n \varphi_{n-1}(x)\f$, where \f$b_n = \frac{n^2}{4n^2-1}\f$
 *
 * so the Legendre polynomials are the characteristic polynomial of the tridiagonal Jacobi matrix
 *
 * \f$J = \begin{pmatrix}
 * 0 & \sqrt{b_1} \\
 * \sqrt{b_1} & 0 & \sqrt{b_2} \\
 * & \ddots &  & \ddots \\
 *  & & \sqrt{b_{N-2}} & 0 & \sqrt{b_{N-1}} \\
 *  & & & \sqrt{b_{N-1}} & 0 \end{pmatrix}\f$
 *
 * and therefore the quadrature locations are the eigenvalues of \f$J\f$
 */
template <class scalar_t>
Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> gauss_impl(size_t N)
{
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> mat_t;
	mat_t J(N, N);
	J.setZero();

	// Fill the diagonals of the matrix
	for (size_t n = 0; n < N-1; ++n)
	{
		scalar_t v = scalar_t(n+1);
		scalar_t u = v / sqrt(4.0*v*v - 1);
		J(n,n+1) = u;
		J(n+1,n) = u;
	}

	// Get the eigenvalues and eigenvectors
	Eigen::SelfAdjointEigenSolver<mat_t> solver(J);
	Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> V;
	V.resize(N,2);
	V.col(0) = solver.eigenvalues();
	V.col(1) = 2.0 * solver.eigenvectors().row(0).cwiseAbs2();

	return V;
}

// forward declaration
template <class Domain>
class gaussian_quadrature;


/**
 * \brief traits of a Gaussian quadrature
 * \tparam Domain the quadrature domain
 */
template <class Domain>
struct quadrature_traits<gaussian_quadrature<Domain> >
{
	/** \brief type of the domain */
	typedef Domain domain_t;
};

/**
 * \brief Gaussian quadrature over a line domain
 */
template <>
class gaussian_quadrature<line_domain> :
	public quadrature_base<gaussian_quadrature<line_domain> >
{
public:
	/** \brief the base class */
	typedef quadrature_base<gaussian_quadrature<line_domain> > base_t;
	/** \brief the location type */
	typedef base_t::xi_t xi_t;
	/** \brief the scalar type*/
	typedef base_t::scalar_t scalar_t;

	/** \brief self-returning */
	typedef gaussian_quadrature type;

	/**
	 * \brief default constructor creating an empty quadrature
	 */
	gaussian_quadrature() :
		base_t(0)
	{
	}

	/**
	 * \brief constructor for a given polynomial degree
	 * \param [in] degree polynomial degree of the quadrature
	 */
	gaussian_quadrature(size_t degree) :
		base_t(degree/2+1)
	{
		size_t N = degree/2+1;
		// compute 1D Gaussian locations and weights
		auto V = gauss_impl<scalar_t>(N);

		// Fill the points and weights
		for(size_t i = 0; i < N; ++i)
		{
			xi_t xi;
			xi << V(i,0);
			push_back(quadrature_elem_t(xi, V(i,1)));
		}
	}
};


/**
 * \brief Gaussian quadrature over a quad domain
 */
template <>
class gaussian_quadrature<quad_domain> :
	public quadrature_base<gaussian_quadrature<quad_domain> >
{
public:
	/** \brief the base class */
	typedef quadrature_base<gaussian_quadrature<quad_domain> > base_t;
	/** \brief the domain type */
	typedef base_t::domain_t domain_t;
	/** \brief the location type */
	typedef base_t::xi_t xi_t;
	/** \brief the scalar type */
	typedef base_t::scalar_t scalar_t;

	/** \brief self-returning */
	typedef gaussian_quadrature type;

	/**
	 * \brief default constructor creating an empty quadrature
	 */
	gaussian_quadrature() :
		base_t(0)
	{
	}

	/**
	 * \brief constructor for a given polynomial order
	 * \param degree polynomial order
	 */
	gaussian_quadrature(size_t degree) :
		base_t((degree/2+1) * (degree/2+1))
	{
		size_t N = degree/2+1;
		auto V = gauss_impl<scalar_t>(N);

		// Fill the points and weights
		for(size_t i = 0; i < N; ++i)
		{
			for(size_t j = 0; j < N; ++j)
			{
				xi_t xi;
				xi << V(i,0), V(j,0);
				push_back(quadrature_elem_t(xi, V(i,1)*V(j,1)));
			}
		}
	}
}; // end of class gauss_quad


/**
 * \brief number of quadrature points for different Dunavant orders
 */
static size_t const dunavant_num[] = {
	1 /*0*/,
	1 /*1*/,
	3 /*2*/,
	4 /*3*/,
	6 /*4*/,
	7 /*5*/,
	12 /*6*/,
	13 /*7*/,
	16 /*8*/,
	19 /*9*/,
	23 /*10*/,
	33 /*11*/,
	33 /*12*/,
	37 /*13*/,
	42 /*14*/
	};

/**
 * \brief specialisation of gauss_quadrature for a triangle domain
 */
template <>
class gaussian_quadrature<tria_domain> :
	public quadrature_base<gaussian_quadrature<tria_domain> >
{
public:
	/** \brief base class */
	typedef quadrature_base<gaussian_quadrature<tria_domain> > base_t;
	/** \brief the quadrature elem */
	typedef base_t::quadrature_elem_t quadrature_elem_t;
	/** \brief the quadrature location type */
	typedef base_t::xi_t xi_t;

	/** \brief self-returning */
	typedef gaussian_quadrature type;


	/**
	 * \brief default constructor creating an empty quadrature
	 */
	gaussian_quadrature() :
		base_t(0)
	{
	}

	/**
	 * \brief constructor for a given polynomial order
	 * \param degree polynomial order
	 */
	gaussian_quadrature(size_t degree) :
		base_t(dunavant_num[degree])
	{
		switch(degree)
		{
		case 0: case 1:
			push_back(quadrature_elem_t(xi_t(1./3., 1./3.0), 1./2.0));
			break;
		case 2:
			push_back(quadrature_elem_t(xi_t(1./6., 4./6.0), 1./6.0));
			push_back(quadrature_elem_t(xi_t(1./6., 1./6.0), 1./6.0));
			push_back(quadrature_elem_t(xi_t(4./6., 1./6.0), 1./6.0));
			break;
		case 3:
			push_back(quadrature_elem_t(xi_t(1./3., 1./3.0), -0.281250000000000));
			push_back(quadrature_elem_t(xi_t(1./5., 3./5.0),  0.260416666666667));
			push_back(quadrature_elem_t(xi_t(1./5., 1./5.0),  0.260416666666667));
			push_back(quadrature_elem_t(xi_t(3./5., 1./5.0),  0.260416666666667));
			break;
		case 4:
			push_back(quadrature_elem_t(xi_t(0.445948490915965, 0.108103018168070),  0.111690794839006));
			push_back(quadrature_elem_t(xi_t(0.445948490915965, 0.445948490915965),  0.111690794839006));
			push_back(quadrature_elem_t(xi_t(0.108103018168070, 0.445948490915965),  0.111690794839006));
			push_back(quadrature_elem_t(xi_t(0.091576213509771, 0.816847572980459),  0.054975871827661));
			push_back(quadrature_elem_t(xi_t(0.091576213509771, 0.091576213509771),  0.054975871827661));
			push_back(quadrature_elem_t(xi_t(0.816847572980459, 0.091576213509771),  0.054975871827661));
			break;
		case 5:
			push_back(quadrature_elem_t(xi_t(0.333333333333333, 0.333333333333333),  0.112500000000000));
			push_back(quadrature_elem_t(xi_t(0.470142064105115, 0.059715871789770),  0.066197076394253));
			push_back(quadrature_elem_t(xi_t(0.470142064105115, 0.470142064105115),  0.066197076394253));
			push_back(quadrature_elem_t(xi_t(0.059715871789770, 0.470142064105115),  0.066197076394253));
			push_back(quadrature_elem_t(xi_t(0.101286507323456, 0.797426985353087),  0.062969590272414));
			push_back(quadrature_elem_t(xi_t(0.101286507323456, 0.101286507323456),  0.062969590272414));
			push_back(quadrature_elem_t(xi_t(0.797426985353087, 0.101286507323456),  0.062969590272414));
			break;
		case 6:
			push_back(quadrature_elem_t(xi_t(0.249286745170910, 0.501426509658179),  0.058393137863189));
			push_back(quadrature_elem_t(xi_t(0.249286745170910, 0.249286745170910),  0.058393137863189));
			push_back(quadrature_elem_t(xi_t(0.501426509658179, 0.249286745170910),  0.058393137863189));
			push_back(quadrature_elem_t(xi_t(0.063089014491502, 0.873821971016996),  0.025422453185104));
			push_back(quadrature_elem_t(xi_t(0.063089014491502, 0.063089014491502),  0.025422453185104));
			push_back(quadrature_elem_t(xi_t(0.873821971016996, 0.063089014491502),  0.025422453185104));
			push_back(quadrature_elem_t(xi_t(0.310352451033784, 0.053145049844817),  0.041425537809187));
			push_back(quadrature_elem_t(xi_t(0.636502499121399, 0.310352451033784),  0.041425537809187));
			push_back(quadrature_elem_t(xi_t(0.053145049844817, 0.636502499121399),  0.041425537809187));
			push_back(quadrature_elem_t(xi_t(0.053145049844817, 0.310352451033784),  0.041425537809187));
			push_back(quadrature_elem_t(xi_t(0.310352451033784, 0.636502499121399),  0.041425537809187));
			push_back(quadrature_elem_t(xi_t(0.636502499121399, 0.053145049844817),  0.041425537809187));
			break;
		case 7:
			push_back(quadrature_elem_t(xi_t(0.333333333333333, 0.333333333333333), -0.074785022233841));
			push_back(quadrature_elem_t(xi_t(0.260345966079040, 0.479308067841920), 0.087807628716604));
			push_back(quadrature_elem_t(xi_t(0.260345966079040, 0.260345966079040), 0.087807628716604));
			push_back(quadrature_elem_t(xi_t(0.479308067841920, 0.260345966079040), 0.087807628716604));
			push_back(quadrature_elem_t(xi_t(0.065130102902216, 0.869739794195568), 0.026673617804419));
			push_back(quadrature_elem_t(xi_t(0.065130102902216, 0.065130102902216), 0.026673617804419));
			push_back(quadrature_elem_t(xi_t(0.869739794195568, 0.065130102902216), 0.026673617804419));
			push_back(quadrature_elem_t(xi_t(0.312865496004874, 0.048690315425316), 0.038556880445128));
			push_back(quadrature_elem_t(xi_t(0.638444188569810, 0.312865496004874), 0.038556880445128));
			push_back(quadrature_elem_t(xi_t(0.048690315425316, 0.638444188569810), 0.038556880445128));
			push_back(quadrature_elem_t(xi_t(0.048690315425316, 0.312865496004874), 0.038556880445128));
			push_back(quadrature_elem_t(xi_t(0.312865496004874, 0.638444188569810), 0.038556880445128));
			push_back(quadrature_elem_t(xi_t(0.638444188569810, 0.048690315425316), 0.038556880445128));
			break;
		case 8:
			push_back(quadrature_elem_t(xi_t(0.333333333333333 ,  0.333333333333333), 0.072157803838894));
			push_back(quadrature_elem_t(xi_t(0.459292588292723 ,  0.081414823414554), 0.047545817133642));
			push_back(quadrature_elem_t(xi_t(0.459292588292723 ,  0.459292588292723), 0.047545817133642));
			push_back(quadrature_elem_t(xi_t(0.081414823414554 ,  0.459292588292723), 0.047545817133642));
			push_back(quadrature_elem_t(xi_t(0.170569307751760 ,  0.658861384496480), 0.051608685267359));
			push_back(quadrature_elem_t(xi_t(0.170569307751760 ,  0.170569307751760), 0.051608685267359));
			push_back(quadrature_elem_t(xi_t(0.658861384496480 ,  0.170569307751760), 0.051608685267359));
			push_back(quadrature_elem_t(xi_t(0.050547228317031 ,  0.898905543365938), 0.016229248811599));
			push_back(quadrature_elem_t(xi_t(0.050547228317031 ,  0.050547228317031), 0.016229248811599));
			push_back(quadrature_elem_t(xi_t(0.898905543365938 ,  0.050547228317031), 0.016229248811599));
			push_back(quadrature_elem_t(xi_t(0.263112829634638 ,  0.008394777409958), 0.013615157087218));
			push_back(quadrature_elem_t(xi_t(0.728492392955404 ,  0.263112829634638), 0.013615157087218));
			push_back(quadrature_elem_t(xi_t(0.008394777409958 ,  0.728492392955404), 0.013615157087218));
			push_back(quadrature_elem_t(xi_t(0.008394777409958 ,  0.263112829634638), 0.013615157087218));
			push_back(quadrature_elem_t(xi_t(0.263112829634638 ,  0.728492392955404), 0.013615157087218));
			push_back(quadrature_elem_t(xi_t(0.728492392955404 ,  0.008394777409958), 0.013615157087218));
			break;
		case 9:
			push_back(quadrature_elem_t(xi_t(0.333333333333333,   0.333333333333333), 0.048567898141400));
			push_back(quadrature_elem_t(xi_t(0.489682519198738,   0.020634961602525), 0.015667350113570));
			push_back(quadrature_elem_t(xi_t(0.489682519198738,   0.489682519198738), 0.015667350113570));
			push_back(quadrature_elem_t(xi_t(0.020634961602525,   0.489682519198738), 0.015667350113570));
			push_back(quadrature_elem_t(xi_t(0.437089591492937,   0.125820817014127), 0.038913770502387));
			push_back(quadrature_elem_t(xi_t(0.437089591492937,   0.437089591492937), 0.038913770502387));
			push_back(quadrature_elem_t(xi_t(0.125820817014127,   0.437089591492937), 0.038913770502387));
			push_back(quadrature_elem_t(xi_t(0.188203535619033,   0.623592928761935), 0.039823869463605));
			push_back(quadrature_elem_t(xi_t(0.188203535619033,   0.188203535619033), 0.039823869463605));
			push_back(quadrature_elem_t(xi_t(0.623592928761935,   0.188203535619033), 0.039823869463605));
			push_back(quadrature_elem_t(xi_t(0.044729513394453,   0.910540973211095), 0.012788837829349));
			push_back(quadrature_elem_t(xi_t(0.044729513394453,   0.044729513394453), 0.012788837829349));
			push_back(quadrature_elem_t(xi_t(0.910540973211095,   0.044729513394453), 0.012788837829349));
			push_back(quadrature_elem_t(xi_t(0.221962989160766,   0.036838412054736), 0.021641769688645));
			push_back(quadrature_elem_t(xi_t(0.741198598784498,   0.221962989160766), 0.021641769688645));
			push_back(quadrature_elem_t(xi_t(0.036838412054736,   0.741198598784498), 0.021641769688645));
			push_back(quadrature_elem_t(xi_t(0.036838412054736,   0.221962989160766), 0.021641769688645));
			push_back(quadrature_elem_t(xi_t(0.221962989160766,   0.741198598784498), 0.021641769688645));
			push_back(quadrature_elem_t(xi_t(0.741198598784498,   0.036838412054736), 0.021641769688645));
			break;
		case 10:
			push_back(quadrature_elem_t(xi_t(0.333333333333333,   0.333333333333333), 0.045408995191377));
			push_back(quadrature_elem_t(xi_t(0.028844733232685,   0.485577633383657), 0.018362978878233));
			push_back(quadrature_elem_t(xi_t(0.485577633383657,   0.485577633383657), 0.018362978878233));
			push_back(quadrature_elem_t(xi_t(0.485577633383657,   0.028844733232685), 0.018362978878233));
			push_back(quadrature_elem_t(xi_t(0.781036849029926,   0.109481575485037), 0.022660529717764));
			push_back(quadrature_elem_t(xi_t(0.109481575485037,   0.109481575485037), 0.022660529717764));
			push_back(quadrature_elem_t(xi_t(0.109481575485037,   0.781036849029926), 0.022660529717764));
			push_back(quadrature_elem_t(xi_t(0.141707219414880,   0.307939838764121), 0.036378958422710));
			push_back(quadrature_elem_t(xi_t(0.307939838764121,   0.550352941820999), 0.036378958422710));
			push_back(quadrature_elem_t(xi_t(0.550352941820999,   0.141707219414880), 0.036378958422710));
			push_back(quadrature_elem_t(xi_t(0.307939838764121,   0.141707219414880), 0.036378958422710));
			push_back(quadrature_elem_t(xi_t(0.550352941820999,   0.307939838764121), 0.036378958422710));
			push_back(quadrature_elem_t(xi_t(0.141707219414880,   0.550352941820999), 0.036378958422710));
			push_back(quadrature_elem_t(xi_t(0.025003534762686,   0.246672560639903), 0.014163621265528));
			push_back(quadrature_elem_t(xi_t(0.246672560639903,   0.728323904597411), 0.014163621265528));
			push_back(quadrature_elem_t(xi_t(0.728323904597411,   0.025003534762686), 0.014163621265528));
			push_back(quadrature_elem_t(xi_t(0.246672560639903,   0.025003534762686), 0.014163621265528));
			push_back(quadrature_elem_t(xi_t(0.728323904597411,   0.246672560639903), 0.014163621265528));
			push_back(quadrature_elem_t(xi_t(0.025003534762686,   0.728323904597411), 0.014163621265528));
			push_back(quadrature_elem_t(xi_t(0.009540815400299,   0.066803251012200), 0.004710833481867));
			push_back(quadrature_elem_t(xi_t(0.066803251012200,   0.923655933587500), 0.004710833481867));
			push_back(quadrature_elem_t(xi_t(0.923655933587500,   0.009540815400299), 0.004710833481867));
			push_back(quadrature_elem_t(xi_t(0.066803251012200,   0.009540815400299), 0.004710833481867));
			push_back(quadrature_elem_t(xi_t(0.923655933587500,   0.066803251012200), 0.004710833481867));
			push_back(quadrature_elem_t(xi_t(0.009540815400299,   0.923655933587500), 0.004710833481867));
			break;
		case 11: 	// NOTE: 11 same as 12 because of outliers at order 11
		case 12:
			push_back(quadrature_elem_t(xi_t(0.023565220452390,   0.488217389773805), 0.012865533220227));
			push_back(quadrature_elem_t(xi_t(0.488217389773805,   0.488217389773805), 0.012865533220227));
			push_back(quadrature_elem_t(xi_t(0.488217389773805,   0.023565220452390), 0.012865533220227));
			push_back(quadrature_elem_t(xi_t(0.120551215411079,   0.439724392294460), 0.021846272269019));
			push_back(quadrature_elem_t(xi_t(0.439724392294460,   0.439724392294460), 0.021846272269019));
			push_back(quadrature_elem_t(xi_t(0.439724392294460,   0.120551215411079), 0.021846272269019));
			push_back(quadrature_elem_t(xi_t(0.457579229975768,   0.271210385012116), 0.031429112108943));
			push_back(quadrature_elem_t(xi_t(0.271210385012116,   0.271210385012116), 0.031429112108943));
			push_back(quadrature_elem_t(xi_t(0.271210385012116,   0.457579229975768), 0.031429112108943));
			push_back(quadrature_elem_t(xi_t(0.744847708916828,   0.127576145541586), 0.017398056465355));
			push_back(quadrature_elem_t(xi_t(0.127576145541586,   0.127576145541586), 0.017398056465355));
			push_back(quadrature_elem_t(xi_t(0.127576145541586,   0.744847708916828), 0.017398056465355));
			push_back(quadrature_elem_t(xi_t(0.957365299093579,   0.021317350453210), 0.003083130525780));
			push_back(quadrature_elem_t(xi_t(0.021317350453210,   0.021317350453210), 0.003083130525780));
			push_back(quadrature_elem_t(xi_t(0.021317350453210,   0.957365299093579), 0.003083130525780));
			push_back(quadrature_elem_t(xi_t(0.115343494534698,   0.275713269685514), 0.020185778883191));
			push_back(quadrature_elem_t(xi_t(0.275713269685514,   0.608943235779788), 0.020185778883191));
			push_back(quadrature_elem_t(xi_t(0.608943235779788,   0.115343494534698), 0.020185778883191));
			push_back(quadrature_elem_t(xi_t(0.275713269685514,   0.115343494534698), 0.020185778883191));
			push_back(quadrature_elem_t(xi_t(0.608943235779788,   0.275713269685514), 0.020185778883191));
			push_back(quadrature_elem_t(xi_t(0.115343494534698,   0.608943235779788), 0.020185778883191));
			push_back(quadrature_elem_t(xi_t(0.022838332222257,   0.281325580989940), 0.011178386601152));
			push_back(quadrature_elem_t(xi_t(0.281325580989940,   0.695836086787803), 0.011178386601152));
			push_back(quadrature_elem_t(xi_t(0.695836086787803,   0.022838332222257), 0.011178386601152));
			push_back(quadrature_elem_t(xi_t(0.281325580989940,   0.022838332222257), 0.011178386601152));
			push_back(quadrature_elem_t(xi_t(0.695836086787803,   0.281325580989940), 0.011178386601152));
			push_back(quadrature_elem_t(xi_t(0.022838332222257,   0.695836086787803), 0.011178386601152));
			push_back(quadrature_elem_t(xi_t(0.025734050548330,   0.116251915907597), 0.008658115554329));
			push_back(quadrature_elem_t(xi_t(0.116251915907597,   0.858014033544073), 0.008658115554329));
			push_back(quadrature_elem_t(xi_t(0.858014033544073,   0.025734050548330), 0.008658115554329));
			push_back(quadrature_elem_t(xi_t(0.116251915907597,   0.025734050548330), 0.008658115554329));
			push_back(quadrature_elem_t(xi_t(0.858014033544073,   0.116251915907597), 0.008658115554329));
			push_back(quadrature_elem_t(xi_t(0.025734050548330,   0.858014033544073), 0.008658115554329));
			break;
		case 13:
			push_back(quadrature_elem_t(xi_t(0.333333333333333,   0.333333333333333), 0.026260461700401));
			push_back(quadrature_elem_t(xi_t(0.009903630120591,   0.495048184939705), 0.005640072604665));
			push_back(quadrature_elem_t(xi_t(0.495048184939705,   0.495048184939705), 0.005640072604665));
			push_back(quadrature_elem_t(xi_t(0.495048184939705,   0.009903630120591), 0.005640072604665));
			push_back(quadrature_elem_t(xi_t(0.062566729780852,   0.468716635109574), 0.015711759181227));
			push_back(quadrature_elem_t(xi_t(0.468716635109574,   0.468716635109574), 0.015711759181227));
			push_back(quadrature_elem_t(xi_t(0.468716635109574,   0.062566729780852), 0.015711759181227));
			push_back(quadrature_elem_t(xi_t(0.170957326397447,   0.414521336801277), 0.023536251252097));
			push_back(quadrature_elem_t(xi_t(0.414521336801277,   0.414521336801277), 0.023536251252097));
			push_back(quadrature_elem_t(xi_t(0.414521336801277,   0.170957326397447), 0.023536251252097));
			push_back(quadrature_elem_t(xi_t(0.541200855914337,   0.229399572042831), 0.023681793268178));
			push_back(quadrature_elem_t(xi_t(0.229399572042831,   0.229399572042831), 0.023681793268178));
			push_back(quadrature_elem_t(xi_t(0.229399572042831,   0.541200855914337), 0.023681793268178));
			push_back(quadrature_elem_t(xi_t(0.771151009607340,   0.114424495196330), 0.015583764522897));
			push_back(quadrature_elem_t(xi_t(0.114424495196330,   0.114424495196330), 0.015583764522897));
			push_back(quadrature_elem_t(xi_t(0.114424495196330,   0.771151009607340), 0.015583764522897));
			push_back(quadrature_elem_t(xi_t(0.950377217273082,   0.024811391363459), 0.003987885732537));
			push_back(quadrature_elem_t(xi_t(0.024811391363459,   0.024811391363459), 0.003987885732537));
			push_back(quadrature_elem_t(xi_t(0.024811391363459,   0.950377217273082), 0.003987885732537));
			push_back(quadrature_elem_t(xi_t(0.094853828379579,   0.268794997058761), 0.018424201364366));
			push_back(quadrature_elem_t(xi_t(0.268794997058761,   0.636351174561660), 0.018424201364366));
			push_back(quadrature_elem_t(xi_t(0.636351174561660,   0.094853828379579), 0.018424201364366));
			push_back(quadrature_elem_t(xi_t(0.268794997058761,   0.094853828379579), 0.018424201364366));
			push_back(quadrature_elem_t(xi_t(0.636351174561660,   0.268794997058761), 0.018424201364366));
			push_back(quadrature_elem_t(xi_t(0.094853828379579,   0.636351174561660), 0.018424201364366));
			push_back(quadrature_elem_t(xi_t(0.018100773278807,   0.291730066734288), 0.008700731651911));
			push_back(quadrature_elem_t(xi_t(0.291730066734288,   0.690169159986905), 0.008700731651911));
			push_back(quadrature_elem_t(xi_t(0.690169159986905,   0.018100773278807), 0.008700731651911));
			push_back(quadrature_elem_t(xi_t(0.291730066734288,   0.018100773278807), 0.008700731651911));
			push_back(quadrature_elem_t(xi_t(0.690169159986905,   0.291730066734288), 0.008700731651911));
			push_back(quadrature_elem_t(xi_t(0.018100773278807,   0.690169159986905), 0.008700731651911));
			push_back(quadrature_elem_t(xi_t(0.022233076674090,   0.126357385491669), 0.007760893419522));
			push_back(quadrature_elem_t(xi_t(0.126357385491669,   0.851409537834241), 0.007760893419522));
			push_back(quadrature_elem_t(xi_t(0.851409537834241,   0.022233076674090), 0.007760893419522));
			push_back(quadrature_elem_t(xi_t(0.126357385491669,   0.022233076674090), 0.007760893419522));
			push_back(quadrature_elem_t(xi_t(0.851409537834241,   0.126357385491669), 0.007760893419522));
			push_back(quadrature_elem_t(xi_t(0.022233076674090,   0.851409537834241), 0.007760893419522));
			break;
		case 14:
			push_back(quadrature_elem_t(xi_t(0.022072179275643,   0.488963910362179), 0.010941790684715));
			push_back(quadrature_elem_t(xi_t(0.488963910362179,   0.488963910362179), 0.010941790684715));
			push_back(quadrature_elem_t(xi_t(0.488963910362179,   0.022072179275643), 0.010941790684715));
			push_back(quadrature_elem_t(xi_t(0.164710561319092,   0.417644719340454), 0.016394176772063));
			push_back(quadrature_elem_t(xi_t(0.417644719340454,   0.417644719340454), 0.016394176772063));
			push_back(quadrature_elem_t(xi_t(0.417644719340454,   0.164710561319092), 0.016394176772063));
			push_back(quadrature_elem_t(xi_t(0.453044943382323,   0.273477528308839), 0.025887052253646));
			push_back(quadrature_elem_t(xi_t(0.273477528308839,   0.273477528308839), 0.025887052253646));
			push_back(quadrature_elem_t(xi_t(0.273477528308839,   0.453044943382323), 0.025887052253646));
			push_back(quadrature_elem_t(xi_t(0.645588935174913,   0.177205532412543), 0.021081294368497));
			push_back(quadrature_elem_t(xi_t(0.177205532412543,   0.177205532412543), 0.021081294368497));
			push_back(quadrature_elem_t(xi_t(0.177205532412543,   0.645588935174913), 0.021081294368497));
			push_back(quadrature_elem_t(xi_t(0.876400233818255,   0.061799883090873), 0.007216849834889));
			push_back(quadrature_elem_t(xi_t(0.061799883090873,   0.061799883090873), 0.007216849834889));
			push_back(quadrature_elem_t(xi_t(0.061799883090873,   0.876400233818255), 0.007216849834889));
			push_back(quadrature_elem_t(xi_t(0.961218077502598,   0.019390961248701), 0.002461701801200));
			push_back(quadrature_elem_t(xi_t(0.019390961248701,   0.019390961248701), 0.002461701801200));
			push_back(quadrature_elem_t(xi_t(0.019390961248701,   0.961218077502598), 0.002461701801200));
			push_back(quadrature_elem_t(xi_t(0.057124757403648,   0.172266687821356), 0.012332876606282));
			push_back(quadrature_elem_t(xi_t(0.172266687821356,   0.770608554774996), 0.012332876606282));
			push_back(quadrature_elem_t(xi_t(0.770608554774996,   0.057124757403648), 0.012332876606282));
			push_back(quadrature_elem_t(xi_t(0.172266687821356,   0.057124757403648), 0.012332876606282));
			push_back(quadrature_elem_t(xi_t(0.770608554774996,   0.172266687821356), 0.012332876606282));
			push_back(quadrature_elem_t(xi_t(0.057124757403648,   0.770608554774996), 0.012332876606282));
			push_back(quadrature_elem_t(xi_t(0.092916249356972,   0.336861459796345), 0.019285755393531));
			push_back(quadrature_elem_t(xi_t(0.336861459796345,   0.570222290846683), 0.019285755393531));
			push_back(quadrature_elem_t(xi_t(0.570222290846683,   0.092916249356972), 0.019285755393531));
			push_back(quadrature_elem_t(xi_t(0.336861459796345,   0.092916249356972), 0.019285755393531));
			push_back(quadrature_elem_t(xi_t(0.570222290846683,   0.336861459796345), 0.019285755393531));
			push_back(quadrature_elem_t(xi_t(0.092916249356972,   0.570222290846683), 0.019285755393531));
			push_back(quadrature_elem_t(xi_t(0.014646950055654,   0.298372882136258), 0.007218154056767));
			push_back(quadrature_elem_t(xi_t(0.298372882136258,   0.686980167808088), 0.007218154056767));
			push_back(quadrature_elem_t(xi_t(0.686980167808088,   0.014646950055654), 0.007218154056767));
			push_back(quadrature_elem_t(xi_t(0.298372882136258,   0.014646950055654), 0.007218154056767));
			push_back(quadrature_elem_t(xi_t(0.686980167808088,   0.298372882136258), 0.007218154056767));
			push_back(quadrature_elem_t(xi_t(0.014646950055654,   0.686980167808088), 0.007218154056767));
			push_back(quadrature_elem_t(xi_t(0.001268330932872,   0.118974497696957), 0.002505114419250));
			push_back(quadrature_elem_t(xi_t(0.118974497696957,   0.879757171370171), 0.002505114419250));
			push_back(quadrature_elem_t(xi_t(0.879757171370171,   0.001268330932872), 0.002505114419250));
			push_back(quadrature_elem_t(xi_t(0.118974497696957,   0.001268330932872), 0.002505114419250));
			push_back(quadrature_elem_t(xi_t(0.879757171370171,   0.118974497696957), 0.002505114419250));
			push_back(quadrature_elem_t(xi_t(0.001268330932872,   0.879757171370171), 0.002505114419250));
			break;
		default:
			throw std::out_of_range("unsupported dunavant degree");
			break;
		}
	}
};

/**
 * \brief tag for the family of Gaussian quadratures
 */
struct gauss_family_tag { typedef gauss_family_tag type; };

/** \brief specialisation of quadrature_type to Gaussian family on line */
template<class Domain>
struct quadrature_type<gauss_family_tag, Domain> :
	gaussian_quadrature<Domain>
{
};

}

#endif // GAUSSIAN_QUADRATURE_HPP_INCLUDED
