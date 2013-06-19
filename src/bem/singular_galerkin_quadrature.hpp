/**
 * \file singular_galerkin_quadrature.hpp
 * \brief implementation of singular Galerkin quadratures
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */
#ifndef SINGULAR_GALERKIN_QUADRATURE_HPP_INCLUDED
#define SINGULAR_GALERKIN_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"

/** \brief singularity types */
enum singularity_type {
	REGULAR,		/**< \brief no singularity */
	FACE_MATCH,	/**< \brief two elements are identical */
	EDGE_MATCH,	/**< \brief two elements share common edge */
	CORNER_MATCH	/**< \brief two elements share common corner */
};


/**
 * \brief class computing singular Galerkin type quadratures for different domains
 * \tparam quadrature_family_t the regular quadrature family
 * \tparam test_domain_t the test domain type
 * \tparam trial_domain_t the trial domain type
 */
template <class quadrature_family_t, class test_domain_t, class trial_domain_t>
class singular_galerkin_quadrature;



/**
 * \brief base structure for quadrature helpers
 * \tparam test_domain_t type of the test domain
 * \tparam trial_domain_t type of the trial domain
 */
template <class test_domain_t, class trial_domain_t>
struct helper_base
{
	/** \brief the scalar type of the singular quadrature */
	typedef typename test_domain_t::scalar_t scalar_t;
	/** \brief the matrix type of a Descartes quadrature element */
	typedef Eigen::Matrix<scalar_t,
		test_domain_t::dimension + trial_domain_t::dimension, 1> descartes_quad_t;
};

/**
 * \brief helper structure for the tria-tria case
 * \tparam match_type the type of singulartiy
 */
template <singularity_type match_type>
class tria_helper;

/**
 * \brief specialisation of ::tria_helper for the FACE_MATCH case
 */
template <>
class tria_helper<FACE_MATCH>
  : public helper_base<tria_domain, tria_domain>
{
public:
	/** \brief indicates whether the quadrature described below is symmetric or not */
	static bool const is_symmetric = true;
	/** \brief the number of subdomains of the present singular quadrature */
	static unsigned const num_domains = 3;
	
	/**
	 * \brief transform the 4d quadrature into a singular Duffy one
	 * \param [in,out] x the quadrature points
	 * \param [in,out] w the quadrature weight
	 * \param [in] idx the subdomain id
	 */
	static void transform(descartes_quad_t &x, scalar_t &w, unsigned idx)
	{
		scalar_t mu1, mu2, xi1, xi2, J;
		switch (idx)
		{
		case 0:
			mu1 = x(0);
			mu2 = x(0) * x(1);
			xi1 = (1-mu1) * x(2);
			xi2 = xi1 * x(3);
			J = (1-mu1) * xi1;
			break;
		case 1:
			mu1 = x(0) * x(1);
			mu2 = x(0) * (x(1) - 1.0);
			xi1 = (1-mu1+mu2) * x(2) - mu2;
			xi2 = (xi1+mu2) * x(3) - mu2;
			J = (1.0-mu1+mu2) * (xi1+mu2);
			break;
		case 2:
			mu1 = x(0) * x(1);
			mu2 = x(0);
			xi1 = (1-mu2) * x(2) + mu2 - mu1;
			xi2 = (xi1-mu2+mu1) * x(3);
			J = (1-mu2) * (xi1-mu2+mu1);
			break;
		}

		w *= J * x(0);

		// transform the quadratures back to the (0,0) (1,0) (0,1) simplex
		x(0) = -xi1 + 1.0;
		x(1) = xi2;
		x(2) = -(mu1+xi1)+1.0;
		x(3) = mu2+xi2;
	}
};

/**
 * \brief specialisation of ::tria_helper for the EDGE_MATCH case
 */
template <>
class tria_helper<EDGE_MATCH>
  : public helper_base<tria_domain, tria_domain>
{
public:
	/** \brief indicates whether the quadrature described below is symmetric or not */
	static bool const is_symmetric = false;
	/** \brief the number of subdomains of the present singular quadrature */
	static unsigned const num_domains = 6;

	/**
	 * \brief transform the 4d quadrature into a singular Duffy one
	 * \param [in,out] x the quadrature points
	 * \param [in,out] w the quadrature weight
	 * \param [in] idx the subdomain id
	 */
	static void transform(descartes_quad_t &x, scalar_t &w, unsigned idx)
	{
		scalar_t mu1, mu2, xi1, xi2;
		switch (idx)
		{
		case 0:
			mu1 = -x(0) * x(1);
			mu2 = -x(0) * x(1) * x(2);
			xi1 = (1.0-x(0)) * x(3) + x(0);
			xi2 = x(0) * (1.0-x(1) + x(1)*x(2));
			break;
		case 1:
			mu1 = x(0) * x(1);
			mu2 = x(0) * x(1) * x(2);
			xi1 = (1.0-x(0)) * x(3) + x(0) * (1.0-x(1));
			xi2 = x(0) * (1.0-x(1));
			break;
		case 2:
			mu1 = -x(0) * x(1) * x(2);
			mu2 = x(0) * x(1) * (1.0 - x(2));
			xi1 = (1.0-x(0)) * x(3) + x(0);
			xi2 = x(0) * (1.0-x(1));
			break;
		case 3:
			mu1 = x(0) * x(1) * x(2);
			mu2 = x(0) * x(1) * (x(2) - 1.0);
			xi1 = (1.0-x(0)) * x(3) + x(0) * (1.0 - x(1)*x(2));
			xi2 = x(0) * (1.0 - x(1)*x(2));
			break;
		case 4:
			mu1 = -x(0) * x(1) * x(2);
			mu2 = -x(0) * x(1);
			xi1 = (1.0-x(0)) * x(3) + x(0);
			xi2 = x(0);
			break;
		case 5:
			mu1 = x(0) * x(1) * x(2);
			mu2 = x(0) * x(1);
			xi1 = (1.0-x(0)) * x(3) + x(0) * (1.0-x(1)*x(2));
			xi2 = x(0) * (1.0-x(1));
			break;
		}
		double J = x(1) * x(0)*x(0) * (1.0-x(0));
		w *= J;

		// transform the quadratures back to the (0,0) (1,0) (0,1) simplex
		x(0) = -xi1 + 1.0;
		x(1) = xi2;
		x(2) = -(mu1+xi1)+1.0;
		x(3) = mu2+xi2;
	}
};

/**
 * \brief specialisation of ::tria_helper for the CORNER_MATCH case
 */
template <>
class tria_helper<CORNER_MATCH>
  : public helper_base<tria_domain, tria_domain>
{
public:
	/** \brief indicates whether the quadrature described below is symmetric or not */
	static bool const is_symmetric = true;
	/** \brief the number of subdomains of the present singular quadrature */
	static unsigned const num_domains = 1;
	
	/**
	 * \brief transform the 4d quadrature into a singular Duffy one
	 * \param [in,out] x the quadrature points
	 * \param [in,out] w the quadrature weight
	 */
	static void transform(descartes_quad_t &x, scalar_t &w, unsigned)
	{
		scalar_t xi1 = x(0);
		scalar_t xi2 = xi1*x(1);
		scalar_t eta1 = xi1*x(2);
		scalar_t eta2 = eta1*x(3);
		scalar_t J = x(2)*x(0)*x(0)*x(0);

		w *= J;

		// transform the quadratures back to the (0,0) (1,0) (0,1) simplex
		x(0) = -xi1 + 1.0;
		x(1) = xi2;
		x(2) = -eta1+1.0;
		x(3) = eta2;
	}
};


/**
 * \brief specialisation of ::singular_galerkin_quadrature for the tria-tria case
 * \details The implementation follows the paper of Tattaglia and Barzini
 * \tparam quadrature_family_t the regular quadrature family
 */
template <class quadrature_family_t>
class singular_galerkin_quadrature<quadrature_family_t, tria_domain, tria_domain>
{
public:
	/** \brief the (regular) quadrature type */
	typedef typename quadrature_type<quadrature_family_t, tria_domain>::type quadrature_t;
	/** \brief the quadrature element type */
	typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;

	/**
	 * \brief generate a singular quadrature for a given singularity type
	 * \tparam match_type the singularity type
	 * \param [out] test_quadrature the test quadrature to be extended
	 * \param [out] trial_quadrature the trial quadrature to be extended
	 * \param [in] SINGULARITY_ORDER polynomial order of the underlying regular quadrature
	 * \todo traversing four line quadratures should be replaced by traversing their
	 * Descartes product
	 */
	template <singularity_type match_type>
	static void generate(
		quadrature_t &test_quadrature,
		quadrature_t &trial_quadrature,
		unsigned SINGULARITY_ORDER)
	{
		/** \brief the helper class type */
		typedef tria_helper<match_type> hlp_t;
		/** \brief the scalar type */
		typedef typename hlp_t::scalar_t scalar_t;
		/** \brief the double quadrature element type */
		typedef typename hlp_t::descartes_quad_t descartes_quad_t;

		/** \brief the regular line quadrature type */
		typename quadrature_type<quadrature_family_t, line_domain>::type base_quad(SINGULARITY_ORDER);
		// transform the regular line quadrature into the (0,1) domain
		Eigen::Matrix<scalar_t, 2, 1> c(0.0, 1.0);
		base_quad.template transform_inplace<line_1_shape_set>(c);

		// traversing the four dimensional regular quadrature elements
		for (unsigned i1 = 0; i1 < base_quad.size(); ++i1)
		{
			scalar_t x1 = base_quad[i1].get_xi()(0);
			scalar_t w1 = base_quad[i1].get_w();
			for (unsigned i2 = 0; i2 < base_quad.size(); ++i2)
			{
				scalar_t x2 = base_quad[i2].get_xi()(0);
				scalar_t w2 = base_quad[i2].get_w();
				for (unsigned i3 = 0; i3 < base_quad.size(); ++i3)
				{
					scalar_t x3 = base_quad[i3].get_xi()(0);
					scalar_t w3 = base_quad[i3].get_w();
					for (unsigned i4 = 0; i4 < base_quad.size(); ++i4)
					{
						scalar_t x4 = base_quad[i4].get_xi()(0);
						scalar_t w4 = base_quad[i4].get_w();

						for (unsigned idx = 0; idx < hlp_t::num_domains; ++idx)
						{
							descartes_quad_t x(x1, x2, x3, x4);
							scalar_t w = w1 * w2 * w3 * w4;
							
							// transform the 4d quadrature into the desired one
							hlp_t::transform(x, w, idx);
							
							// separate into test and trial quadratures
							test_quadrature.push_back(quadrature_elem_t(x.topRows(2), w));
							trial_quadrature.push_back(quadrature_elem_t(x.bottomRows(2), 1.0));

							// and vica versa if symmetry requires
							if (hlp_t::is_symmetric)
							{
								test_quadrature.push_back(quadrature_elem_t(x.bottomRows(2), w));
								trial_quadrature.push_back(quadrature_elem_t(x.topRows(2), 1.0));
							}
						}
					} // for loop for i4
				} // for loop for i3
			} // for loop for i2
		} // for loop for i1
	} // function generate
};


template <singularity_type match_type>
struct quad_helper;

template <>
struct quad_helper<FACE_MATCH> : helper_base<quad_domain, quad_domain>
{
	static const unsigned num_domains = 8;
	static const bool is_symmetric = false;
	
	int corners[num_domains][4][2] = {
		{{0, 0}, {0, 0}, { 1,  0}, { 1,  1}},
		{{0, 0}, {0, 0}, { 1,  1}, { 0,  1}},
		{{0, 0}, {0, 0}, { 0,  1}, {-1,  1}},
		{{0, 0}, {0, 0}, {-1,  1}, {-1,  0}},
		{{0, 0}, {0, 0}, {-1,  0}, {-1, -1}},
		{{0, 0}, {0, 0}, {-1, -1}, { 0, -1}},
		{{0, 0}, {0, 0}, { 0, -1}, { 1, -1}},
		{{0, 0}, {0, 0}, { 1, -1}, { 1,  0}}
	};
};


/**
 * \brief specialisation of ::singular_galerkin_quadrature for the quad-quad case
 * \details The implementation follows our algorithm
 * \tparam quadrature_family_t the regular quadrature family
 */
template <class quadrature_family_t>
class singular_galerkin_quadrature<quadrature_family_t, quad_domain, quad_domain>
{
public:
	/** \brief the (regular) quadrature type */
	typedef typename quadrature_type<quadrature_family_t, quad_domain>::type quadrature_t;
	/** \brief the quadrature element type */
	typedef typename quadrature_t::quadrature_elem_t quadrature_elem_t;

	/**
	 * \brief generate a singular quadrature for a given singularity type
	 * \tparam match_type the singularity type
	 * \param [out] test_quadrature the test quadrature to be extended
	 * \param [out] trial_quadrature the trial quadrature to be extended
	 * \param [in] SINGULARITY_ORDER polynomial order of the underlying regular quadrature
	 */
	template <singularity_type match_type>
	static void generate(
		quadrature_t &test_quadrature,
		quadrature_t &trial_quadrature,
		unsigned SINGULARITY_ORDER)
	{
		/** \brief the helper class type */
		typedef quad_helper<match_type> hlp_t;

		// create a regular quad quadrature for Duffy transformation purposes
		quadrature_t base_quad(SINGULARITY_ORDER);
		
		test_quadrature = base_quad;
		trial_quadrature = base_quad;
	}
};



#endif
