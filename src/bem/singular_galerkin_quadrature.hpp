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
	static const unsigned num_domains = 4;
	static const bool is_symmetric = true;
	static quad_domain::scalar_t const corners[4][4][2];

	typedef quad_domain::xi_t xi_t;
	static xi_t transform_eta(xi_t const &eta)
	{
		return eta;
	}
};

quad_domain::scalar_t
	const quad_helper<FACE_MATCH>::corners[4][4][2] = {
		{{0, 0}, {0, 0}, { 2, -2}, { 2,  0}},
		{{0, 0}, {0, 0}, { 2,  0}, { 2,  2}},
		{{0, 0}, {0, 0}, { 2,  2}, { 0,  2}},
		{{0, 0}, {0, 0}, { 0,  2}, {-2,  2}}
};

template <>
struct quad_helper<EDGE_MATCH> : helper_base<quad_domain, quad_domain>
{
	static const unsigned num_domains = 6;
	static const bool is_symmetric = false;
	static quad_domain::scalar_t const corners[6][4][2];

	typedef quad_domain::xi_t xi_t;
	static xi_t transform_eta(xi_t const &eta)
	{
		return xi_t(eta(0), -2.0-eta(1));
	}
};

quad_domain::scalar_t
	const quad_helper<EDGE_MATCH>::corners[6][4][2] = {
		{{ 0,  0}, { 0,  0}, {-2,  0}, {-2, -2}},
		{{ 0,  0}, { 0,  0}, {-2, -2}, { 0, -2}},
		{{ 0,  0}, { 0,  0}, { 0, -2}, { 2, -2}},
		{{ 0,  0}, { 0,  0}, { 2, -2}, { 2,  0}},
		{{-2, -2}, {-2, -4}, { 0, -4}, { 0, -2}},
		{{ 0, -2}, { 0, -4}, { 2, -4}, { 2, -2}},
};

template <>
struct quad_helper<CORNER_MATCH> : helper_base<quad_domain, quad_domain>
{
	static const unsigned num_domains = 5;
	static const bool is_symmetric = false;
	static quad_domain::scalar_t const corners[5][4][2];

	typedef quad_domain::xi_t xi_t;
	static xi_t transform_eta(xi_t const &eta)
	{
		return -2.0 * xi_t::Ones() - eta;
	}
};

quad_domain::scalar_t
	const quad_helper<CORNER_MATCH>::corners[5][4][2] = {
		{{ 0,  0}, { 0,  0}, {-2,  0}, {-2, -2}},
		{{ 0,  0}, { 0,  0}, {-2, -2}, { 0, -2}},
		{{-4,  0}, {-4, -2}, {-2, -2}, {-2,  0}},
		{{-4, -2}, {-4, -4}, {-2, -4}, {-2, -2}},
		{{-2, -2}, {-2, -4}, { 0, -4}, { 0, -2}}
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
		typedef typename hlp_t::scalar_t scalar_t;

		// create a regular quad quadrature for Duffy transformation purposes
		quadrature_t base_quad(SINGULARITY_ORDER);

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
				eta_lims[0] = hlp_t::transform_eta(quad_domain::get_corners()[0]);
				eta_lims[1] = hlp_t::transform_eta(quad_domain::get_corners()[2]);

				// make sure that the Y rectangle has positive side lengths
				for (unsigned j = 0; j < 2; ++j)
					if (eta_lims[0](j) > eta_lims[1](j))
						std::swap(eta_lims[0](j), eta_lims[1](j));

				// compute (-mu + Y) * X
				xi_t mu = outer_quad[out_idx].get_xi();
				scalar_t out_w = outer_quad[out_idx].get_w();
				corners <<
					-mu(0)+eta_lims[0](0), -mu(1)+eta_lims[0](1),
					-mu(0)+eta_lims[1](0), -mu(1)+eta_lims[0](1),
					-mu(0)+eta_lims[1](0), -mu(1)+eta_lims[1](1),
					-mu(0)+eta_lims[0](0), -mu(1)+eta_lims[1](1);
				for (int i = 0; i < corners.rows(); ++i)
					for (int j = 0; j < corners.cols(); ++j)
						/** \todo -1.0 and +1.0 come from the domain corners, this hard coding is sick */
							corners(i,j) = std::max(std::min(corners(i,j), 1.0), -1.0);

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
	/** \brief the (regular) quadrature type */
	typedef typename quadrature_type<quadrature_family_t, quad_domain>::type test_quadrature_t;
	typedef typename quadrature_type<quadrature_family_t, tria_domain>::type trial_quadrature_t;
	/** \brief the quadrature element type */
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;
	/** \brief location type of the outer and inner quadratures */
	typedef typename quadrature_elem_t::xi_t xi_t;

	typedef singular_galerkin_quadrature<quadrature_family_t, tria_domain, tria_domain> base_sing_t;

	/**
	* \brief generate a singular quadrature for a given singularity type
	* \tparam match_type the singularity type
	* \param [out] test_quadrature the test quadrature to be extended
	* \param [out] trial_quadrature the trial quadrature to be extended
	* \param [in] SINGULARITY_ORDER polynomial order of the underlying regular quadrature
	*/
	template <singularity_type match_type>
	static void generate(
		test_quadrature_t &test_quadrature,
		trial_quadrature_t &trial_quadrature,
		unsigned SINGULARITY_ORDER);

	template <>
	static void generate<CORNER_MATCH>(
		test_quadrature_t &test_quadrature,
		trial_quadrature_t &trial_quadrature,
		unsigned SINGULARITY_ORDER)
	{
		trial_quadrature_t test_base;
		trial_quadrature_t trial_base;
		base_sing_t::template generate<CORNER_MATCH>(
			test_base, trial_base, SINGULARITY_ORDER);

		int idx[2][3] = {
			{0, 1, 2},
			{0, 2, 3}
		};

		for (unsigned d = 0; d < 2; ++d)
		{
			Eigen::Matrix<tria_domain::scalar_t, 3, 2> corners;
			for (unsigned i = 0; i < 3; ++i)
				corners.row(i) = quad_domain::get_corners()[idx[d][i]].transpose();
			trial_quadrature_t test_trans = test_base.template transform<tria_1_shape_set>(corners);

			for (unsigned i = 0; i < test_trans.size(); ++i)
			{
				test_quadrature.push_back(test_trans[i]);
				trial_quadrature.push_back(trial_base[i]);
			}
		}
	}



	template <>
	static void generate<EDGE_MATCH>(
		test_quadrature_t &test_quadrature,
		trial_quadrature_t &trial_quadrature,
		unsigned SINGULARITY_ORDER)
	{
		Eigen::Matrix<tria_domain::scalar_t, 3, 2> corners;
		trial_quadrature_t test_base;
		trial_quadrature_t trial_base;

		base_sing_t::template generate<EDGE_MATCH>(
			test_base, trial_base, SINGULARITY_ORDER);

		corners << -1.0, -1.0, 1.0, -1.0, 1.0, 1.0;
		test_base.template transform_inplace<tria_1_shape_set>(corners);
		for (unsigned i = 0; i < test_base.size(); ++i)
		{
			test_quadrature.push_back(test_base[i]);
			trial_quadrature.push_back(trial_base[i]);
		}

		test_base.clear();
		trial_base.clear();

		singular_galerkin_quadrature<quadrature_family_t, tria_domain, tria_domain>::template generate<CORNER_MATCH>(
			test_base, trial_base, SINGULARITY_ORDER);

		corners << -1.0, -1.0, 1.0, 1.0, -1.0, 1.0;
		test_base.template transform_inplace<tria_1_shape_set>(corners);
		for (unsigned i = 0; i < test_base.size(); ++i)
		{
			test_quadrature.push_back(test_base[i]);
			trial_quadrature.push_back(trial_base[i]);
		}
	}
};


template <class quadrature_family_t>
class singular_galerkin_quadrature<quadrature_family_t, tria_domain, quad_domain>
{
public:
	typedef typename quadrature_type<quadrature_family_t, tria_domain>::type test_quadrature_t;
	typedef typename quadrature_type<quadrature_family_t, quad_domain>::type trial_quadrature_t;

	/**
	* \brief generate a singular quadrature for a given singularity type
	* \tparam match_type the singularity type
	* \param [out] test_quadrature the test quadrature to be extended
	* \param [out] trial_quadrature the trial quadrature to be extended
	* \param [in] SINGULARITY_ORDER polynomial order of the underlying regular quadrature
	*/
	template <singularity_type match_type>
	static void generate(
		test_quadrature_t &test_quadrature,
		trial_quadrature_t &trial_quadrature,
		unsigned SINGULARITY_ORDER)
	{
		singular_galerkin_quadrature<quadrature_family_t, quad_domain, tria_domain>::generate<match_type>(
			trial_quadrature, test_quadrature, SINGULARITY_ORDER);
	}
};


#endif
