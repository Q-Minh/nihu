#ifndef SINGULAR_GALERKIN_QUADRATURE_HPP_INCLUDED
#define SINGULAR_GALERKIN_QUADRATURE_HPP_INCLUDED

#include "quadrature.hpp"

/** \brief singularity types */
enum singularity_type {
	REGULAR,	/**< \brief no singularity */
	FACE_MATCH,	/**< \brief two elements are identical */
	EDGE_MATCH,	/**< \brief two elements share common edge */
	CORNER_MATCH	/**< \brief two elements share common corner */
};


template <class quadrature_family_t, class test_domain_t, class trial_domain_t>
class singular_galerkin_quadrature
{
public:
	typedef typename quadrature_type<quadrature_family_t, test_domain_t>::type test_quadrature_t;
	typedef typename quadrature_type<quadrature_family_t, trial_domain_t>::type trial_quadrature_t;
	typedef typename test_quadrature_t::quadrature_elem_t quadrature_elem_t;
	typedef typename quadrature_elem_t::scalar_t scalar_t;
	typedef Eigen::Matrix<scalar_t, test_domain_t::dimension+trial_domain_t::dimension, 1> quadr4_t;

	template <singularity_type match_type, class test_dom_t, class trial_dom_t>
	class helper;

	template <>
	class helper<FACE_MATCH, tria_domain, tria_domain>
	{
	public:
		static bool const is_symmetric = true;
		static unsigned const num_domains = 3;

		static void transform(quadr4_t &x, scalar_t &w, unsigned idx)
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

	template <>
	class helper<EDGE_MATCH, tria_domain, tria_domain>
	{
	public:
		static bool const is_symmetric = false;
		static unsigned const num_domains = 6;

		static void transform(quadr4_t &x, scalar_t &w, unsigned idx)
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

	template <>
	class helper<CORNER_MATCH, tria_domain, tria_domain>
	{
	public:
		static bool const is_symmetric = true;
		static unsigned const num_domains = 1;

		static void transform(quadr4_t &x, scalar_t &w, unsigned idx)
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

	template <singularity_type match_type>
	static void generate(
		test_quadrature_t &m_face_test_quadrature,
		trial_quadrature_t &m_face_trial_quadrature,
		unsigned SINGULARITY_ORDER)
	{
		typedef typename helper<match_type, test_domain_t, trial_domain_t> inter;

		typename quadrature_type<quadrature_family_t, line_domain>::type base_quad(SINGULARITY_ORDER);
		Eigen::Matrix<scalar_t, 2, 1> c(0.0, 1.0);
		base_quad.transform_inplace<line_1_shape_set>(c);

		for (unsigned i1 = 0; i1 < base_quad.size(); ++i1)
		{
			scalar_t x1 = base_quad[i1].get_xi()(0,0);
			scalar_t w1 = base_quad[i1].get_w();
			for (unsigned i2 = 0; i2 < base_quad.size(); ++i2)
			{
				scalar_t x2 = base_quad[i2].get_xi()(0,0);
				scalar_t w2 = base_quad[i2].get_w();
				for (unsigned i3 = 0; i3 < base_quad.size(); ++i3)
				{
					scalar_t x3 = base_quad[i3].get_xi()(0,0);
					scalar_t w3 = base_quad[i3].get_w();
					for (unsigned i4 = 0; i4 < base_quad.size(); ++i4)
					{
						scalar_t x4 = base_quad[i4].get_xi()(0,0);
						scalar_t w4 = base_quad[i4].get_w();

						for (unsigned idx = 0; idx < inter::num_domains; ++idx)
						{
							quadr4_t x = quadr4_t(x1, x2, x3, x4);
							scalar_t w = w1 * w2 * w3 * w4;
							inter::transform(x, w, idx);

							m_face_test_quadrature.push_back(quadrature_elem_t(x.topRows(2), w));
							m_face_trial_quadrature.push_back(quadrature_elem_t(x.bottomRows(2), 1.0));

							if (inter::is_symmetric)
							{
								m_face_test_quadrature.push_back(quadrature_elem_t(x.bottomRows(2), w));
								m_face_trial_quadrature.push_back(quadrature_elem_t(x.topRows(2), 1.0));
							}
						}
					} // for loop for i4
				} // for loop for i3
			} // for loop for i2
		} // for loop for i1
	} // function generate
};

#endif
