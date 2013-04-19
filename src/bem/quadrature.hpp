/**
 * \file quadrature.hpp
 * \author Peter Fiala fiala@hit.bme.hu Peter Rucz rucz@hit.bme.hu
 * \brief implementation of class ::quadrature_elem, ::quadrature_base and gaussian quadratures
 */
#ifndef QUADRATURE_HPP_INCLUDED
#define QUADRATURE_HPP_INCLUDED

#include "shapeset.hpp"

#include <type_traits> // is_same

#include <iostream>
#include <stdexcept>

#include <Eigen/StdVector>
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >


/**
 * \brief a quadrature element is a base point and a weight
 * \tparam XiType the local coordinate type
 * \tparam ScalarType the scalar of given coordinates
 */
template <class XiType, class ScalarType>
class quadrature_elem
{
public:
	typedef XiType xi_t;			/**< \brief template argument as tested type */
	typedef ScalarType scalar_t;	/**< \brief template argument as nested type */

	/**
	 * \brief constructor initialising all members
	 * \param [in] xi base location
	 * \param [in] w weight
	 */
	quadrature_elem(xi_t const &xi = xi_t(), scalar_t const &w = scalar_t())
		: m_xi(xi), m_w(w)
	{
	}

	/**
	 * \brief return constant reference to base point
	 * \return reference to base point
	 */
	xi_t const &get_xi(void) const
	{
		return m_xi;
	}

	/**
	 * \brief return constant reference to weight
	 * \return reference to weight
	 */
	scalar_t const &get_w(void) const
	{
		return m_w;
	}

	/**
	 * \brief transform a location according to a shape function set
	 * \tparam LSet the shape function set
	 * \param [in] coords transformation corners
	 * \return reference to the transformed object
	 */
	template <class LSet>
	quadrature_elem &transform_inplace(const Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> &coords)
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// compute Jacobian and multiply with original weight
		m_w *= (LSet::eval_dshape(m_xi).transpose() * coords).determinant();
		// compute new quadrature location
		m_xi = LSet::eval_shape(m_xi).transpose() * coords;
		return *this;
	}

	// struct is used in a std::vector, therefore this declaration is necessary
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
	xi_t m_xi;		/**< \brief base point */
	scalar_t m_w;	/**< \brief weight */
};


// forward declaration of singular traits
template <class quadrature>
struct singular_traits;

// forward declaration of quadrature traits
template <class Derived>
struct quadrature_traits;


/**
 * \brief CRTP base class of all quadratures
 * \tparam Derived the CRTP derived class
 */
template <class Derived>
class quadrature_base :
	public std::vector<quadrature_elem<typename quadrature_traits<Derived>::domain_t::xi_t, typename quadrature_traits<Derived>::domain_t::scalar_t> , Eigen::aligned_allocator<quadrature_elem<typename quadrature_traits<Derived>::domain_t::xi_t, typename quadrature_traits<Derived>::domain_t::scalar_t> > >
{
protected:
	/**
	 * \brief static cast to derived const type
	 */
	Derived const &derived(void) const
	{
		return static_cast<Derived const &>(*this);
	}

	/**
	 * \brief static cast to derived type
	 */
	Derived &derived(void)
	{
		return static_cast<Derived &>(*this);
	}

public:
	typedef quadrature_traits<Derived> traits_t;				/**< \brief corresponding traits type */
	typedef typename traits_t::domain_t domain_t;				/**< \brief the domain type */
	typedef typename domain_t::xi_t xi_t; 						/**< \brief local coordinate type */
	typedef typename domain_t::scalar_t scalar_t; 				/**< \brief local scalar type */
	typedef quadrature_elem<xi_t, scalar_t> quadrature_elem_t;	/**< \brief the quadrature elem type */

	/**
	 * \brief constructor allocating space for the quadrature elements
	 * \param N number of quadrature elements
	 */
	quadrature_base(unsigned N = 0)
	{
		this->reserve(N);
	}

	/**
	 * \brief print a quadrature into an output stream
	 * \param os the output stream
	 * \return reference to the stream
	 */
	std::ostream & print (std::ostream & os) const
	{
		os << "Base points: \t Weights:" << std::endl;
		std::for_each(this->begin(), this->end(), [&os] (quadrature_elem_t const &e) {
			os << e.get_xi().transpose() << "\t\t" << e.get_w() << std::endl;
		});
		os << "Sum of weights: " <<
			std::accumulate(this->begin(), this->end(), scalar_t(), [] (scalar_t x, quadrature_elem_t const &qe) {
			return x + qe.get_w();
		}) << std::endl;

		return os;
	}

	/**
	 * \brief transform the domain of the quadrature with a given shape set and corner points
	 * \tparam LSet shape set type of transformation
	 * \param [in] coords corner coordinates of transformed domain
	 * \return transformed quadrature
	 */
	template <class LSet>
	Derived transform(const Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> &coords) const
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// dimension check
		static_assert(std::is_same<xi_t, typename LSet::xi_t>::value,
			"Quadrature and shape set dimensions must match");
		Derived result(derived());
		return result.transform_inplace<LSet>(coords);
	}

	/**
	 * \brief transform the domain of the quadrature in place
	 * \tparam LSet shape set type of transformation
	 * \param [in] coords corner coordinates of transformed domain
	 * \return reference to the transformed quadrature
	 */
	template <class LSet>
	Derived &transform_inplace(const Eigen::Matrix<scalar_t, LSet::num_nodes, LSet::domain_t::dimension> &coords)
	{
		// CRTP check
		static_assert(std::is_base_of<shape_set_base<LSet>, LSet>::value,
			"LSet must be derived from shape_set_base<LSet>");
		// dimension check
		static_assert(std::is_same<xi_t, typename LSet::xi_t>::value,
			"Quadrature and shape set dimensions must match");
		for (auto it = this->begin(); it != this->end(); ++it)
			it->transform_inplace<LSet>(coords);
		return derived();
	}

	/**
	 * \brief add two quadratures
	 * \details It is assured that the dimensions match.
	 * \tparam otherDerived other quadrature type
	 * \param [in] other the other quadrature
	 * \return new quadrature
	 */
	template <class otherDerived>
	Derived operator +(const quadrature_base<otherDerived> &other) const
	{
		static_assert(std::is_same<xi_t, typename otherDerived::xi_t>::value,
			"Quadrature domain dimensions must match");
		Derived result(this->size()+other.size());
		result.insert(result.end(), this->begin(), this->end());
		result.insert(result.end(), other.begin(), other.end());
		return result;
	}

	/**
	 * \brief add another quadrature to this
	 * \details It is assured that the dimensions match.
	 * \tparam otherDerived other quadrature type
	 * \param [in] other the other quadrature
	 * \return reference to the new quadrature
	 */
	template <class otherDerived>
	Derived &operator +=(const quadrature_base<otherDerived> &other)
	{
		static_assert(std::is_same<xi_t, typename otherDerived::xi_t>::value,
			"Quadrature domain dimensions must match");
		this->reserve(this->size()+other.size());
		this->insert(this->end(), other.begin(), other.end());
		return derived();
	}

	/**
	 * \brief create a singular quadrature from a selected internal point
	 * \param [in] degree quadrature degree of the generating element
	 * \param [in] singular_point the internal singular point
	 */
	static Derived singular_quadrature_inside(unsigned degree, xi_t const &singular_point)
	{
		typedef singular_traits<Derived> singular_traits_t;
		typedef typename singular_traits_t::singular_source_type singular_source_t;
		typedef typename singular_traits_t::transformation_lset transformation_lset;
		const unsigned nResolution = singular_traits_t::nResolution;
		const unsigned nCorners = singular_traits_t::nCorners;

		xi_t const *corners = domain_t::get_corners();

		Derived result;
		singular_source_t source_quad(degree);
		for (unsigned r = 0; r < nResolution; ++r)
		{
			Eigen::Matrix<scalar_t, nCorners, domain_t::dimension> coords;
			for (unsigned c = 0; c < nCorners; ++c)
			{
				if (singular_traits_t::index_inside[r][c] == -1)
					coords.row(c) = singular_point;
				else
					coords.row(c) = corners[singular_traits_t::index_inside[r][c]];
			}
			result += source_quad.transform<transformation_lset>(coords);
		}
		return result;
	}
}; // quadrature_base


/**
 * \brief print a quadrature into an ouput stream
 * \tparam the quadrature's domain type
 * \param os the output stream
 * \param Q the quadrature
 * \return the modified output stream
 */
template<class Derived>
std::ostream & operator << (std::ostream & os, const quadrature_base<Derived>& Q)
{
	return Q.print(os);
}

/**
 * \brief return 1D N-point Guassian quadrature
 * \tparam scalar_t the scalar type
 * \param N number of points
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
 * \brief Gaussian quadrature
 */
class gauss_line;

/**
 * \brief traits of gaussian line quadrature
 */
template <>
struct quadrature_traits<gauss_line>
{
	typedef line_domain domain_t;	/**< \brief type of the domain */
};

/**
 * \brief specialisation of gauss_quadrature for a line domain
 */
class gauss_line : public quadrature_base<gauss_line>
{
public:
	typedef quadrature_base<gauss_line> base_t;	/**< \brief the base class */
	typedef base_t::xi_t xi_t;	/**< \brief the locatin type */
	typedef base_t::scalar_t scalar_t;	/**< \brief the scalar type*/


	/**
	 * \brief default constructor creating an empty quadrature
	 */
	gauss_line() : base_t(0) {};

	/**
	 * \brief constructor for a given polynomial degree
	 * \param degree polynomial degree
	 */
	gauss_line(unsigned degree) : base_t(degree/2+1)
	{
		unsigned N = degree/2+1;
		typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> eig_t;
		eig_t V = gauss_impl<scalar_t>(N);

		// Fill the points and weights
		for(unsigned i = 0; i < N; ++i)
		{
			xi_t xi;
			xi << V(i,0);
			push_back(quadrature_elem_t(xi, V(i,1)));
		}
	}
};


// forward declaration
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
 * \brief specialisation of gauss_quadrature for a quad domain
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
	gauss_quad() : base_t(0) {};

	/**
	 * \brief constructor for a given polynomial order
	 * \param degree polynomial order
	 */
	gauss_quad(unsigned degree) : base_t((degree/2+1) * (degree/2+1))
	{
		unsigned N = degree/2+1;
		typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> eig_t;
		eig_t V = gauss_impl<scalar_t>(N);

		// Fill the points and weights
		for(unsigned i = 0; i < N; ++i)
			for(unsigned j = 0; j < N; ++j)
			{
				xi_t xi;
				xi << V(i,0), V(j,0);
				push_back(quadrature_elem_t(xi, V(i,1)*V(j,1)));
			}
	}
};


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
// Quadrature families

/**
 * \brief tag for the family of gaussian quadratures
 */
struct gauss_family_tag;

/**
 * \brief assign a quadrature type to a quadrature family and a domain
 */
template <class Family, class Domain>
struct quadrature_type;

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

#endif

