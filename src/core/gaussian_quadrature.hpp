/**
 * \file gaussian_quadrature.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief implementation of Gaussian quadratures
 */

#ifndef GAUSSIAN_QUADRATURE_HPP_INCLUDED
#define GAUSSIAN_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"

#include <stdexcept>

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
Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> gauss_impl(unsigned N)
{
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> mat_t;
	mat_t J(N, N);
	J.setZero();

	// Fill the diagonals of the matrix
	for (unsigned n = 0; n < N-1; ++n)
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
	gaussian_quadrature(unsigned degree) :
		base_t(degree/2+1)
	{
		unsigned N = degree/2+1;
		// compute 1D Gaussian locations and weights
		auto V = gauss_impl<scalar_t>(N);

		// Fill the points and weights
		for(unsigned i = 0; i < N; ++i)
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
	gaussian_quadrature(unsigned degree) :
		base_t((degree/2+1) * (degree/2+1))
	{
		unsigned N = degree/2+1;
		auto V = gauss_impl<scalar_t>(N);

		// Fill the points and weights
		for(unsigned i = 0; i < N; ++i)
		{
			for(unsigned j = 0; j < N; ++j)
			{
				xi_t xi;
				xi << V(i,0), V(j,0);
				push_back(quadrature_elem_t(xi, V(i,1)*V(j,1)));
			}
		}
	}
}; // end of class gauss_quad


/**
 * \brief number of quadrature points for different dunavat orders
 */
static unsigned const dunavant_num[] = {1, 1, 3, 4, 6, 7, 12, 13, 16, 19};

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
	gaussian_quadrature(unsigned degree) :
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
		default:
			throw std::out_of_range("unsupported dunavant degree");
			break;
		}
	}
};

/**
 * \brief tag for the family of Gaussian quadratures
 */
struct gauss_family_tag;

/** \brief specialisation of quadrature_type to Gaussian family on line */
template<class Domain>
struct quadrature_type<gauss_family_tag, Domain> :
	gaussian_quadrature<Domain>
{
};

#endif // GAUSSIAN_QUADRATURE_HPP_INCLUDED
