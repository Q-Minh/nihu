// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2013  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2013  Peter Rucz <rucz@hit.bme.hu>
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
 * \file singular_galerkin_quadrature.hpp
 * \ingroup quadrature
 * \brief implementation of singular Galerkin quadratures
 * \details This file implements class ::singular_galerkin_quadrature that generates
 * quadratures integrating weakly singular (1/r) kernels in the Galerkin BEM context.
 */
#ifndef SINGULAR_GALERKIN_QUADRATURE_HPP_INCLUDED
#define SINGULAR_GALERKIN_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"
#include "element_match.hpp"

/**
* \brief class computing singular Galerkin type quadratures for different domains
* \tparam quadrature_family_t the regular quadrature family
* \tparam test_domain_t the test domain type
* \tparam trial_domain_t the trial domain type
* \details The class generates singular quadratures that integrate weakly singular (1/r)
* kernels in the Galerkin BEM context. The class is directly specialised for different
* domain pairs. Each specialisation is based on a slightly different algorithm, but the basic
* approach is the same:
* - Domain variables are transformed into coordinate differences;
* - The resulting complex shaped domains are subdivided into simple subdomains
* - Subdomains are transformed into polar coordinates using a Duffy-type coordinate transform.
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
template <class match_type>
struct tria_helper;

/**
* \brief specialisation of ::tria_helper for the FACE_MATCH case
*/
template <>
struct tria_helper<match::face_match_type>
	: public helper_base<tria_domain, tria_domain>
{
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
	static void transform_inplace(descartes_quad_t &x, scalar_t &w, unsigned idx)
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

		x(0) = xi1;
		x(1) = xi2;
		x(2) = mu1+xi1;	// eta1
		x(3) = mu2+xi2;	// eta2
	}
};

/**
* \brief specialisation of ::tria_helper for the EDGE_MATCH case
*/
template <>
struct tria_helper<match::edge_match_type>
	: public helper_base<tria_domain, tria_domain>
{
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
	static void transform_inplace(descartes_quad_t &x, scalar_t &w, unsigned idx)
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

		x(0) = xi1;
		x(1) = xi2;
		x(2) = mu1+xi1;	// eta1
		x(3) = mu2+xi2;	// eta2
	}
};

/**
* \brief specialisation of ::tria_helper for the CORNER_MATCH case
*/
template <>
struct tria_helper<match::corner_match_type>
	: public helper_base<tria_domain, tria_domain>
{
	/** \brief indicates whether the quadrature described below is symmetric or not */
	static bool const is_symmetric = true;
	/** \brief the number of subdomains of the present singular quadrature */
	static unsigned const num_domains = 1;

	/**
	* \brief transform the 4d quadrature into a singular Duffy one
	* \param [in,out] x the quadrature points
	* \param [in,out] w the quadrature weight
	*/
	static void transform_inplace(descartes_quad_t &x, scalar_t &w, unsigned)
	{
		scalar_t xi1 = x(0);
		scalar_t xi2 = xi1*x(1);
		scalar_t eta1 = xi1*x(2);
		scalar_t eta2 = eta1*x(3);
		scalar_t J = x(2)*x(0)*x(0)*x(0);

		w *= J;

		x(0) = xi1;
		x(1) = xi2;
		x(2) = eta1;
		x(3) = eta2;
	}
};


/**
* \brief specialisation of ::singular_galerkin_quadrature for the tria-tria case
* \tparam quadrature_family_t the regular quadrature family
* \details The implementation follows the paper of Tattaglia and Barzini
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
	* \param [in] singular_quadrature_order polynomial order of the underlying regular quadrature
	* \todo traversing four line quadratures should be replaced by traversing their
	* Descartes product
	* \todo singularity type should use singularity type and template instantiation from parameter
	*/
	template <class match_type>
	static void generate(
		quadrature_t &test_quadrature,
		quadrature_t &trial_quadrature,
		unsigned singular_quadrature_order)
	{
		/** \brief the helper class type */
		typedef tria_helper<match_type> hlp_t;
		/** \brief the scalar type */
		typedef typename hlp_t::scalar_t scalar_t;
		/** \brief the double quadrature element type */
		typedef typename hlp_t::descartes_quad_t descartes_quad_t;

		/** \brief the regular line quadrature type */
		typename quadrature_type<quadrature_family_t, line_domain>::type base_quad(singular_quadrature_order);
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
							// generate the descartes quadrature element
							descartes_quad_t x(x1, x2, x3, x4);
							scalar_t w = w1 * w2 * w3 * w4;

							// transform the 4d quadrature into the desired singular one
							hlp_t::transform_inplace(x, w, idx);

							// transform back into our standard triangle simplex
							x(0) -= x(1);
							x(2) -= x(3); // Jacobian is 1.0

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


/** \brief helper struct of the quad-quad algorithm
 * \tparam match_type the type of singularity
 */
template <class match_type>
struct quad_helper;

/** \brief specialisation of ::quad_helper for the ::FACE_MATCH case */
template <>
struct quad_helper<match::face_match_type> : helper_base<quad_domain, quad_domain>
{
	/** \brief number of subdomains */
	static const unsigned num_domains = 4;
	/** \brief indicates whether quadrature is symmetric to two variables */
	static const bool is_symmetric = true;
	/** \brief type of the domain variable */
	typedef quad_domain::xi_t xi_t;
	/** \brief type of the domain scalar */
	typedef quad_domain::scalar_t scalar_t;
	/** \brief corners of subdomains */
	static scalar_t const corners[4][4][2];
	/** \brief transform vector from Y' domain to Y domain */
	static xi_t transform_eta(xi_t const &eta)
	{
		return eta;
	}
};

quad_domain::scalar_t
	const quad_helper<match::face_match_type>::corners[4][4][2] = {
		{{0.0, 0.0}, {0.0, 0.0}, { 2.0, -2.0}, { 2.0,  0.0}},
		{{0.0, 0.0}, {0.0, 0.0}, { 2.0,  0.0}, { 2.0,  2.0}},
		{{0.0, 0.0}, {0.0, 0.0}, { 2.0,  2.0}, { 0.0,  2.0}},
		{{0.0, 0.0}, {0.0, 0.0}, { 0.0,  2.0}, {-2.0,  2.0}}
};

/** \brief specialisation of ::quad_helper for the ::EDGE_MATCH case */
template <>
struct quad_helper<match::edge_match_type> : helper_base<quad_domain, quad_domain>
{
	/** \brief number of subdomains */
	static const unsigned num_domains = 6;
	/** \brief indicates whether quadrature is symmetric to two variables */
	static const bool is_symmetric = false;
	/** \brief type of the domain variable */
	typedef quad_domain::xi_t xi_t;
	/** \brief type of the domain scalar */
	typedef quad_domain::scalar_t scalar_t;
	/** \brief corners of subdomains */
	static scalar_t const corners[6][4][2];
	/** \brief transform vector from Y' domain to Y domain */
	static xi_t transform_eta(xi_t const &eta)
	{
		return xi_t(eta(0), -2.0-eta(1));
	}
};

quad_domain::scalar_t
	const quad_helper<match::edge_match_type>::corners[6][4][2] = {
		{{ 0.0,  0.0}, { 0.0,  0.0}, {-2.0,  0.0}, {-2.0, -2.0}},
		{{ 0.0,  0.0}, { 0.0,  0.0}, {-2.0, -2.0}, { 0.0, -2.0}},
		{{ 0.0,  0.0}, { 0.0,  0.0}, { 0.0, -2.0}, { 2.0, -2.0}},
		{{ 0.0,  0.0}, { 0.0,  0.0}, { 2.0, -2.0}, { 2.0,  0.0}},
		{{-2.0, -2.0}, {-2.0, -4.0}, { 0.0, -4.0}, { 0.0, -2.0}},
		{{ 0.0, -2.0}, { 0.0, -4.0}, { 2.0, -4.0}, { 2.0, -2.0}},
};

/** \brief specialisation of ::quad_helper for the ::CORNER_MATCH case */
template <>
struct quad_helper<match::corner_match_type> : helper_base<quad_domain, quad_domain>
{
	/** \brief number of subdomains */
	static const unsigned num_domains = 5;
	/** \brief indicates whether quadrature is symmetric to two variables */
	static const bool is_symmetric = false;
	/** \brief type of the domain scalar */
	typedef quad_domain::scalar_t scalar_t;
	/** \brief corners of subdomains */
	static scalar_t const corners[5][4][2];
	/** \brief type of the domain variable */
	typedef quad_domain::xi_t xi_t;
	/** \brief transform vector from Y' domain to Y domain */
	static xi_t transform_eta(xi_t const &eta)
	{
		return -2.0 * xi_t::Ones() - eta;
	}
};

quad_domain::scalar_t
	const quad_helper<match::corner_match_type>::corners[5][4][2] = {
		{{ 0.0,  0.0}, { 0.0,  0.0}, {-2.0,  0.0}, {-2.0, -2.0}},
		{{ 0.0,  0.0}, { 0.0,  0.0}, {-2.0, -2.0}, { 0.0, -2.0}},
		{{-4.0,  0.0}, {-4.0, -2.0}, {-2.0, -2.0}, {-2.0,  0.0}},
		{{-4.0, -2.0}, {-4.0, -4.0}, {-2.0, -4.0}, {-2.0, -2.0}},
		{{-2.0, -2.0}, {-2.0, -4.0}, { 0.0, -4.0}, { 0.0, -2.0}}
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
	/** \brief location type of the outer and inner quadratures */
	typedef typename quadrature_elem_t::xi_t xi_t;

	/**
	* \brief generate a singular quadrature for a given singularity type
	* \tparam match_type the singularity type
	* \param [out] test_quadrature the test quadrature to be extended
	* \param [out] trial_quadrature the trial quadrature to be extended
	* \param [in] singular_quadrature_order polynomial order of the underlying regular quadrature
	*/
	template <class match_type>
	static void generate(
		quadrature_t &test_quadrature,
		quadrature_t &trial_quadrature,
		unsigned singular_quadrature_order)
	{
		/** \brief the helper class type */
		typedef quad_helper<match_type> hlp_t;
		typedef typename hlp_t::scalar_t scalar_t;

		// create a regular quad quadrature for Duffy transformation purposes
		quadrature_t base_quad(singular_quadrature_order);

		for (unsigned d = 0; d < hlp_t::num_domains; ++d)
		{
			// Duffy transformation corners
			Eigen::Matrix<scalar_t, 4, 2> corners;

			// Copy transformation corners from helper struct
			for (unsigned c = 0; c < 4; ++c)
				for (unsigned j = 0; j < 2; ++j)
					corners(c,j) = hlp_t::corners[d][c][j];

			// perform transform
			quadrature_t outer_quad = base_quad.template transform<quad_1_shape_set>(corners);

			for (unsigned out_idx = 0; out_idx < outer_quad.size(); ++out_idx)
			{
				// compute opposite corners of Y rectangle
				xi_t eta_lims[2];
				eta_lims[0] = hlp_t::transform_eta(quad_domain::get_corner(0));
				eta_lims[1] = hlp_t::transform_eta(quad_domain::get_corner(2));

				// make sure that the Y rectangle has positive side lengths
				for (unsigned j = 0; j < 2; ++j)
					if (eta_lims[0](j) > eta_lims[1](j))
						std::swap(eta_lims[0](j), eta_lims[1](j));

				// compute (-mu + Y) intersection with X
				xi_t mu = outer_quad[out_idx].get_xi();
				scalar_t out_w = outer_quad[out_idx].get_w();
				corners <<
					-mu(0)+eta_lims[0](0), -mu(1)+eta_lims[0](1),
					-mu(0)+eta_lims[1](0), -mu(1)+eta_lims[0](1),
					-mu(0)+eta_lims[1](0), -mu(1)+eta_lims[1](1),
					-mu(0)+eta_lims[0](0), -mu(1)+eta_lims[1](1);
				// limit into quad domain
				xi_t const &ximin = quad_domain::get_corner(0);
				xi_t const &ximax = quad_domain::get_corner(2);
				for (int i = 0; i < corners.rows(); ++i)
					for (unsigned j = 0; j < 2; ++j)
						corners(i,j) = std::max(std::min(corners(i,j), ximax(j)), ximin(j));

				quadrature_t inner_quad = base_quad.template transform<quad_1_shape_set>(corners);

				for (unsigned in_idx = 0; in_idx < inner_quad.size(); ++in_idx)
				{
					xi_t xi = inner_quad[in_idx].get_xi();
					xi_t eta = hlp_t::transform_eta(mu + xi);

					scalar_t w = out_w * inner_quad[in_idx].get_w();

					test_quadrature.push_back(quadrature_elem_t(xi, w));
					trial_quadrature.push_back(quadrature_elem_t(eta, 1.0));

					if (hlp_t::is_symmetric)
					{
						test_quadrature.push_back(quadrature_elem_t(eta, 1.0));
						trial_quadrature.push_back(quadrature_elem_t(xi, w));
					}
				} // inner quadrature loop
			} // outer quadrature loop
		} // domain loop
	} // function generate()
};


/**
* \brief specialisation of ::singular_galerkin_quadrature for the quad-tria case
* \details The implementation follows Barzini's algorithm, but the quad is divided into trias
* \tparam quadrature_family_t the regular quadrature family
*/
template <class quadrature_family_t>
class singular_galerkin_quadrature<quadrature_family_t, quad_domain, tria_domain>
{
public:
	/** \brief the (regular) test quadrature type */
	typedef typename quadrature_type<quadrature_family_t, quad_domain>::type test_quadrature_t;
	/** \brief the (regular) trial quadrature type */
	typedef typename quadrature_type<quadrature_family_t, tria_domain>::type trial_quadrature_t;
	/** \brief the quadrature element type */
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;
	/** \brief location type of the outer and inner quadratures */
	typedef typename quadrature_elem_t::xi_t xi_t;
	/** \brief the underlying singular quadrature type (tria-tria) */
	typedef singular_galerkin_quadrature<quadrature_family_t, tria_domain, tria_domain> base_sing_t;

	/**
	* \brief generate a singular quadrature for a given singularity type
	* \tparam match_type the singularity type
	* \param [out] test_quadrature the test quadrature to be extended
	* \param [out] trial_quadrature the trial quadrature to be extended
	* \param [in] singular_quadrature_order polynomial order of the underlying regular quadrature
	*/
	template <class match_type>
	static void generate(
		test_quadrature_t &test_quadrature,
		trial_quadrature_t &trial_quadrature,
		unsigned singular_quadrature_order)
	{
		// call specialised function member
		generate(match_type(), test_quadrature, trial_quadrature, singular_quadrature_order);
	}


private:
	/**
	* \brief specialisation of ::generate for the ::CORNER_MATCH case
	* \param [out] test_quadrature the test quadrature to be extended
	* \param [out] trial_quadrature the trial quadrature to be extended
	* \param [in] singular_quadrature_order polynomial order of the underlying regular quadrature
	*/
	static void generate(
		match::corner_match_type,
		test_quadrature_t &test_quadrature,
		trial_quadrature_t &trial_quadrature,
		unsigned singular_quadrature_order)
	{
		trial_quadrature_t test_base;
		trial_quadrature_t trial_base;
		base_sing_t::template generate<match::corner_match_type>(
			test_base, trial_base, singular_quadrature_order);

		unsigned corner_idx[2][3] = {
			{0, 1, 2},
			{0, 2, 3}
		};

		for (unsigned d = 0; d < 2; ++d)
		{
			Eigen::Matrix<tria_domain::scalar_t, 3, 2> corners;
			for (unsigned i = 0; i < 3; ++i)
				corners.row(i) = quad_domain::get_corner(corner_idx[d][i]);
			trial_quadrature_t test_trans = test_base.template transform<tria_1_shape_set>(corners);

			for (unsigned i = 0; i < test_trans.size(); ++i)
			{
				test_quadrature.push_back(test_trans[i]);
				trial_quadrature.push_back(trial_base[i]);
			}
		}
	}

	/**
	* \brief specialisation of ::generate for the ::EDGE_MATCH case
	* \param [out] test_quadrature the test quadrature to be extended
	* \param [out] trial_quadrature the trial quadrature to be extended
	* \param [in] singular_quadrature_order polynomial order of the underlying regular quadrature
	*/
	static void generate(
		match::edge_match_type,
		test_quadrature_t &test_quadrature,
		trial_quadrature_t &trial_quadrature,
		unsigned singular_quadrature_order)
	{
		// the quad domain is divided into two triangles
		unsigned corner_idx[2][3] = {
			{0, 1, 2},	// EDGE_MATCH
			{0, 2, 3}	// CORNER_MATCH
		};

		// for transformation pusposes
		Eigen::Matrix<tria_domain::scalar_t, 3, 2> corners;

		// the underlying edge-singular tria-tria quadratures
		trial_quadrature_t test_base;
		trial_quadrature_t trial_base;
		base_sing_t::template generate<match::edge_match_type>(
			test_base, trial_base, singular_quadrature_order);

		// assemble corners, transform and insert into result
		for (unsigned i = 0; i < 3; ++i)
			corners.row(i) = quad_domain::get_corner(corner_idx[0][i]);
		test_base.template transform_inplace<tria_1_shape_set>(corners);
		for (unsigned i = 0; i < test_base.size(); ++i)
		{
			test_quadrature.push_back(test_base[i]);
			trial_quadrature.push_back(trial_base[i]);
		}

		// clear the quadrautres before filling them again
		test_base.clear();
		trial_base.clear();

		// the underlying corner-singular tria-tria quadratures
		base_sing_t::template generate<match::corner_match_type>(
			test_base, trial_base, singular_quadrature_order);

		// assemble corners, transform and insert into result
		for (unsigned i = 0; i < 3; ++i)
			corners.row(i) = quad_domain::get_corner(corner_idx[1][i]);
		test_base.template transform_inplace<tria_1_shape_set>(corners);
		for (unsigned i = 0; i < test_base.size(); ++i)
		{
			test_quadrature.push_back(test_base[i]);
			trial_quadrature.push_back(trial_base[i]);
		}
	}	// function generate
};


/**
* \brief specialisation of ::singular_galerkin_quadrature for the tria-quad case
* \details The implementation reuses the quad-tria specialisation
* \tparam quadrature_family_t the regular quadrature family
*/
template <class quadrature_family_t>
class singular_galerkin_quadrature<quadrature_family_t, tria_domain, quad_domain>
{
public:
	/** \brief the regular test quadrature type */
	typedef typename quadrature_type<quadrature_family_t, tria_domain>::type test_quadrature_t;
	/** \brief the regular trial quadrature type */
	typedef typename quadrature_type<quadrature_family_t, quad_domain>::type trial_quadrature_t;

	/**
	* \brief generate a singular quadrature for a given singularity type
	* \tparam match_type the singularity type
	* \param [out] test_quadrature the test quadrature to be extended
	* \param [out] trial_quadrature the trial quadrature to be extended
	* \param [in] singular_quadrature_order polynomial order of the underlying regular quadrature
	*/
	template <class match_type>
	static void generate(
		test_quadrature_t &test_quadrature,
		trial_quadrature_t &trial_quadrature,
		unsigned singular_quadrature_order)
	{
		// call quad-tria version with swapped arguments
		singular_galerkin_quadrature<quadrature_family_t, quad_domain, tria_domain>::
			template generate<match_type>(trial_quadrature, test_quadrature, singular_quadrature_order);
	}
};


#endif

