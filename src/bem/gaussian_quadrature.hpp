/**
 * \file gaussian_quadrature.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief implementation of Gaussian quadratures
 */

 /**
 * \todo Just for fun, an alternative of the Gaussian family should be trapezoid family or something similar.
 */


#ifndef GAUSSIAN_QUADRATURE_HPP_INCLUDED
#define GAUSSIAN_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"

/**
 * \brief return 1D N-point Guassian quadrature
 * \tparam scalar_t the scalar type
 * \param [in] N number of Gaussian points
 * \return matrix containing the Gaussian locations and weights
 */
template <class scalar_t>
Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> gauss_impl(unsigned N)
{
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> mat_t;

	mat_t M(N, N);
	M.setZero();

	// Fill the diagonals of the matrix
	for (unsigned i = 0; i < N-1; ++i)
	{
		scalar_t v = scalar_t(i+1);
		scalar_t u = v / sqrt(4.0*v*v - 1);
		M(i,i+1) = u;
		M(i+1,i) = u;
	}

	// Get the eigenvalues and eigenvectors
	Eigen::SelfAdjointEigenSolver<mat_t> es(M);
	Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> V;
	V.resize(N,2);
	V.col(0) = es.eigenvalues();
	V.col(1) = 2.0 * es.eigenvectors().row(0).cwiseAbs2();

	return V;
}

/**
 * \brief Gaussian quadrature over a line domain
 */
class gauss_line;

/**
 * \brief traits of Gaussian line quadrature
 */
template <>
struct quadrature_traits<gauss_line>
{
	typedef line_domain domain_t;	/**< \brief type of the domain */
};

/**
 * \brief Gaussian quadrature over a line domain
 */
class gauss_line : public quadrature_base<gauss_line>
{
public:
	typedef quadrature_base<gauss_line> base_t;	/**< \brief the base class */
	typedef base_t::xi_t xi_t;	/**< \brief the location type */
	typedef base_t::scalar_t scalar_t;	/**< \brief the scalar type*/

	/**
	 * \brief default constructor creating an empty quadrature
	 */
	gauss_line() : base_t(0)
	{
	}

	/**
	 * \brief constructor for a given polynomial degree
	 * \param [in] degree polynomial degree of the quadrature
	 */
	gauss_line(unsigned degree) : base_t(degree/2+1)
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
}; // class gauss_line


/**
 * \brief Gaussian quadrature over a quad domain
 */
class gauss_quad;

/**
 * \brief traits of gaussian quad quadrature
 */
template <>
struct quadrature_traits<gauss_quad>
{
	typedef quad_domain domain_t;	/**< \brief type of the domain */
};


/**
 * \brief Gaussian quadrature over a quad domain
 */
class gauss_quad : public quadrature_base<gauss_quad>
{
public:
	typedef quadrature_base<gauss_quad> base_t;	/**< \brief the base class */
	typedef base_t::domain_t domain_t;	/**< \brief the domain type */
	typedef base_t::xi_t xi_t;	/**< \brief the location type */
	typedef base_t::scalar_t scalar_t;	/**< \brief the scalar type */

	/**
	 * \brief default constructor creating an empty quadrature
	 */
	gauss_quad() : base_t(0)
	{
	}

	/**
	 * \brief constructor for a given polynomial order
	 * \param degree polynomial order
	 */
	gauss_quad(unsigned degree) : base_t((degree/2+1) * (degree/2+1))
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


// forward declaration
class gauss_tria;


/**
 * \brief traits of gaussian tria quadrature
 */
template <>
struct quadrature_traits<gauss_tria>
{
	typedef tria_domain domain_t;	/**< \brief type of the domain */
};


/**
 * \brief number of quadrature points for different dunavat orders
 */
static unsigned const dunavant_num[] = {1, 1, 3, 4, 6, 7, 12, 13, 16, 19};

/**
 * \brief specialisation of gauss_quadrature for a triangle domain
 */
class gauss_tria : public quadrature_base<gauss_tria>
{
public:
	typedef quadrature_base<gauss_tria> base_t;	/**< \brief base class */
	typedef base_t::quadrature_elem_t quadrature_elem_t;	/**< \brief the quadrature elem */
	typedef base_t::xi_t xi_t; /**< \brief the quadrature location type */


	/**
	 * \brief default constructor creating an empty quadrature
	 */
	gauss_tria() : base_t(0) {}

	/**
	 * \brief constructor for a given polynomial order
	 * \param degree polynomial order
	 */
	gauss_tria(unsigned degree) : base_t(dunavant_num[degree])
	{
		switch(degree)
		{
		case 0:
		case 1:
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


// Quadrature families

/**
 * \brief tag for the family of Gaussian quadratures
 */
struct gauss_family_tag;

/** \brief specialisation of quadrature_type to Gaussian family on line */
template<>
struct quadrature_type<gauss_family_tag, line_domain>
{
	typedef gauss_line type; /**< \brief metafunction return type */
};

/** \brief specialisation of quadrature_type to Gaussian family on triangle */
template<>
struct quadrature_type<gauss_family_tag, tria_domain>
{
	typedef gauss_tria type; /**< \brief metafunction return type */
};

/** \brief specialisation of quadrature_type to Gaussian family on quad */
template<>
struct quadrature_type<gauss_family_tag, quad_domain>
{
	typedef gauss_quad type; /**< \brief metafunction return type */
};



/**
 * \todo the whole singular part from here on should be moved to duffy_quadrature.hpp,
 * and a general singular family type should be introduced that can be specialised to Duffy.
 */


/**
 * \brief singular traits of gaussian quad quadrature
 */
template <>
struct singular_traits<gauss_quad>
{
	/** \brief the singular duffy element is gaussian too */
	typedef gauss_quad singular_source_type;
	/** \brief the tria is divided into duffy quads */
	typedef quad_1_shape_set transformation_lset;
	/** \brief number of corners of a duffy quad = 4 */
	static const unsigned nCorners = singular_source_type::domain_t::num_corners;
	/** \brief number of quads in the division */
	static const unsigned nResolution = 4;
	/** \brief number of quads in the division */
	static const unsigned nCornerResolution = 2;
	/** \brief index matrix describing orientation of the duffy quads */
	static const int index_inside[nResolution][nCorners];
	/** \brief index matrix describing orientation of the duffy quads */
	static const int index_corner[nCornerResolution][nCorners];
};

const int singular_traits<gauss_quad>::index_inside[singular_traits<gauss_quad>::nResolution][singular_traits<gauss_quad>::nCorners] = {
	{0, 1, -1, -1},
	{1, 2, -1, -1},
	{2, 3, -1, -1},
	{3, 0, -1, -1}
};

const int singular_traits<gauss_quad>::index_corner[singular_traits<gauss_quad>::nCornerResolution][singular_traits<gauss_quad>::nCorners] = {
	{0, 0, 1, 2},
	{0, 0, 2, 3}
};

/**
 * \brief singular traits of a gaussian tria quadrature
 */
template <>
struct singular_traits<gauss_tria>
{
	/** \brief the singular duffy element is gaussian too */
	typedef gauss_quad singular_source_type;
	/** \brief duffy quads are transformed with an iso quad L-set */
	typedef quad_1_shape_set transformation_lset;
	/** \brief number of corners of a duffy quad = 4 */
	static const unsigned nCorners = singular_source_type::domain_t::num_corners;
	/** \brief number of triangles in the division */
	static const unsigned nResolution = 3;
	/** \brief number of triangles in the division */
	static const unsigned nCornerResolution = 1;
	/** \brief index matrix describing orientation of the duffy quads */
	static const int index_inside[nResolution][nCorners];
	/** \brief index matrix describing orientation of the duffy quads */
	static const int index_corner[nCornerResolution][nCorners];
};

const int singular_traits<gauss_tria>::index_inside[singular_traits<gauss_tria>::nResolution][singular_traits<gauss_tria>::nCorners] = {
	{0, 1, -1, -1},
	{1, 2, -1, -1},
	{2, 0, -1, -1}
};

const int singular_traits<gauss_tria>::index_corner[singular_traits<gauss_tria>::nCornerResolution][singular_traits<gauss_tria>::nCorners] = {
	{0, 0, 1, 2}
};


#endif // GAUSSIAN_QUADRATURE_HPP_INCLUDED
