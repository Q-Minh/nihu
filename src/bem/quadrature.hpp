#ifndef QUADRATURE_HPP_INCLUDED
#define QUADRATURE_HPP_INCLUDED

#include "domain.hpp"
#include "shapeset.hpp"

#include <iostream>

#include <Eigen/StdVector>
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >

template <class Domain>
class quadrature_elem
{
public:
	typedef Domain domain_t;
	typedef typename domain_t::scalar_t scalar_t;
	typedef typename domain_t::xi_t xi_t;

	quadrature_elem(xi_t const &xi = xi_t(), scalar_t const &w = scalar_t()) : xi(xi), w(w)
	{
	}

	xi_t const &get_xi(void) const
	{
		return xi;
	}

	scalar_t const &get_w(void) const
	{
		return w;
	}

protected:
	xi_t xi;
	scalar_t w;
};

template <class Domain>
class quadrature : public EIGENSTDVECTOR(quadrature_elem<Domain>)
{
public:
	typedef Domain domain_t;
	typedef quadrature_elem<domain_t> quadrature_elem_t;
	typedef typename quadrature_elem_t::scalar_t scalar_t;

	quadrature(unsigned N = 0)
	{
		this->reserve(N);
	}

	std::ostream & print (std::ostream & os) const
	{
		os << "Base points: \t Weights:" << std::endl;
		std::for_each(this->begin(), this->end(), [&os] (quadrature_elem_t const &e) {
			os << e.get_xi() << '\t' << e.get_w() << std::endl;
		});
		return os;
	}

	/*
	template <unsigned nNode>
	quadrature<T, nKer> createSingular(const CShape<T, nKer, nNode>  &LSet, const Matrix<T, nKer, 1> & xi0) const
	{
		quadrature<T, nKer> result(0);

		Matrix<T, nKer, 4> TMatrix;
		const Matrix<T, nKer, nNode> & coords = LSet.getCorners();
		TMatrix.col(0) = TMatrix.col(1) = xi0;

		// Create nNode subquadratures
		for (unsigned k = 0; k < nNode; ++k)
		{
			bool create = true;
			// Fill other nodes
			for (unsigned i = 0; i < 2; ++i)
			{
				TMatrix.col(i+2) = coords.col((i+k)%nNode);
				if ((TMatrix.col(i+2) - xi0).norm() < 0.000001)
				{
					create = false;
					break;
				}
			}
			if (create)
				result+= this->transform<4>(singularTransShape, TMatrix);

		}
		return result;

	}


	// Transformation
	template <unsigned nNode>
	quadrature<T, nKer> transform(const CShape<T, nKer, nNode> &LSet, const Matrix<T, nKer, nNode> &coords) const
	{
		quadrature<T, nKer> result(this->size());
		// Calculate the new local coordinates
		for (unsigned i = 0; i < this->size(); ++i)
		{
			CDescriptor<T, nKer> LD;
			Matrix<T, nKer, 1> xi = (*this)[i].getLocation();
			LD.setLocation(coords*LSet.getShape(xi));
			T jac = abs((coords*LSet.getGradShape(xi)).determinant());
			LD.setWeight((*this)[i].getWeight()*jac);
			result.push_back(LD);
		}
		return result;
	}

	// Transformation in place
	template <unsigned nNode>
	void transformInPlace(const CShape<T, nKer, nNode> &LSet, const Matrix<T, nKer, nNode> &coords)
	{
		for (unsigned i = 0; i < this->size(); ++i)
		{
			// Calculate the new local coordinates
			(*this)[i].setLocation(coords*LSet.getShape((*this)[i].getLocation()));
			// Calculate the jacobian and update weight
			T jac = abs((coords*LSet.getGradShape((*this)[i].getLocation())).determinant());
			(*this)[i].setWeight((*this)[i].getWeight()*jac);
		}
	}

	// Addition
	quadrature<T, nKer> operator +(const quadrature<T, nKer> &other) const
	{
		quadrature<T, nKer> result(size()+other.size());
		result.insert(result.end(), begin(), end());
		result.insert(result.end(), other.begin(), other.end());
		return result;
	}

	quadrature<T, nKer> &operator +=(const quadrature<T, nKer> &other)
	{
		insert(end(), other.begin(), other.end());
		return *this;
	}

	*/
};

/*
template <class T, unsigned nKer>
const CQuad4Shape<T> quadrature<T, nKer>::singularTransShape;
*/

template <class Domain>
class gauss_quadrature;


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


template <>
class gauss_quadrature<line_domain> : public quadrature<line_domain>
{
public:
	typedef quadrature<line_domain> base;
	typedef base::quadrature_elem_t::xi_t xi_t;
	typedef base::scalar_t scalar_t;

	gauss_quadrature(unsigned N) : quadrature<line_domain>(N)
	{
		typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> eig_t;
		eig_t V = gauss_impl<scalar_t>(N);

		// Fill the points and weights
		for(eig_t::Index i = 0; i < N; ++i)
		{
			xi_t xi;
			xi << V(i,0);
			push_back(quadrature_elem_t(xi, V(i,1)));
		}
	}
};


template <>
class gauss_quadrature<quad_domain> : public quadrature<quad_domain>
{
public:
	typedef quadrature<quad_domain> base;
	typedef base::quadrature_elem_t::xi_t xi_t;
	typedef base::scalar_t scalar_t;

	gauss_quadrature(unsigned N) : quadrature<quad_domain>(N*N)
	{
		typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 2> eig_t;
		eig_t V = gauss_impl<scalar_t>(N);

		// Fill the points and weights
		for(eig_t::Index i = 0; i < N; ++i)
		{
			for(eig_t::Index j = 0; j < N; ++j)
			{
				xi_t xi;
				xi << V(i,0), V(j,0);
				push_back(quadrature_elem_t(xi, V(i,1)*V(j,1)));
			}
		}
	}
};

static unsigned const dunavant_num[] = {0, 1, 3, 4, 6, 7, 12, 13, 16, 19};

template<>
class gauss_quadrature<tria_domain> : public quadrature<tria_domain>
{
public:
	typedef quadrature<tria_domain> base;
	typedef base::quadrature_elem_t quadrature_elem_t;
	typedef quadrature_elem_t::xi_t xi_t;

	gauss_quadrature(unsigned degree) : quadrature<tria_domain>(dunavant_num[degree])
	{
		switch(degree)
		{
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
		default:
			break;
		}
	}
};


template<class Domain>
std::ostream & operator << (std::ostream & os, const quadrature<Domain>& Q)
{
	return Q.print(os);
}

#endif
