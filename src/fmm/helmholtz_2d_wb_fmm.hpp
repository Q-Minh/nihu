/** \file helmholtz_2d_wb_fmm.hpp
 * \brief implemenetation of the 2D Wide Band Helmholtz FMM
 */
#ifndef HELMHOLTZ_2D_WB_FMM_HPP_INCLUDED
#define HELMHOLTZ_2D_WB_FMM_HPP_INCLUDED

#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/constants/constants.hpp>

#include "cluster.hpp"
#include "cluster_tree.hpp"
#include "helmholtz_2d_wb_cluster.h"
#include "helmholtz_2d_wb_x2x_matrix.h"
#include "helmholtz_2d_wb_level_data.h"
#include "m2l_indices.hpp"
#include "operator_with_wave_number.hpp"
#include "p2p.hpp"

#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_nearly_singular_integrals.hpp"
#include "library/helmholtz_singular_integrals.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <complex>
#include <iostream>

namespace NiHu
{
namespace fmm
{

/** \brief the 2d wide band helmholetz fmm
 * \tparam WaveNumber the wave number type 
 */
template <class WaveNumber>
class helmholtz_2d_wb_fmm
{
public:
	/** \brief template parameter as nested type */
	typedef WaveNumber wave_number_t;

	typedef helmholtz_2d_wb_cluster cluster_t;
	typedef cluster_tree<cluster_t> cluster_tree_t;
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;
	static size_t const dimension = cluster_t::dimension;
	typedef typename cluster_t::bounding_box_t bounding_box_t;
	typedef typename bounding_box_t::location_t location_t;

	/** \brief the m2m operator */
	class m2m
		: public operator_with_wave_number<wave_number_t>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t cluster_t;
		typedef helmholtz_2d_wb_m2m_matrix result_t;

		m2m(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return bounding_box_t::dist2idx(
				from.get_bounding_box().get_center(),
				to.get_bounding_box().get_center());
		}

		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			using boost::math::cyl_bessel_j;

			location_t const &X = to.get_bounding_box().get_center();
			location_t const &Y = from.get_bounding_box().get_center();
			location_t d = X - Y;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			auto z = this->get_wave_number() * r;
			int Lto = to.get_level_data().get_expansion_length();
			int Lfrom = from.get_level_data().get_expansion_length();
			int L = std::max(Lto, Lfrom);
			cvector_t diag_coeffs(2 * L + 1, 1);
			for (int nu = -L; nu <= L; ++nu)
				diag_coeffs(-nu + L) = std::exp(-std::complex<double>(0., nu*theta))
				* cyl_bessel_j(nu, z);

			return result_t(to.get_level_data(), from.get_level_data(),
				diag_coeffs);
		}
	};

	/** \brief the l2l operator */
	class l2l
		: public operator_with_wave_number<wave_number_t>
	{
	private:
		typedef operator_with_wave_number<WaveNumber> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t cluster_t;
		typedef helmholtz_2d_wb_l2l_matrix result_t;

		l2l(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return bounding_box_t::dist2idx(
				to.get_bounding_box().get_center(),
				from.get_bounding_box().get_center());
		}

		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			using boost::math::cyl_bessel_j;
			using namespace boost::math::double_constants;

			location_t const &X = to.get_bounding_box().get_center();
			location_t const &Y = from.get_bounding_box().get_center();
			location_t d = X - Y;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			auto z = this->get_wave_number() * r;
			int Lto = to.get_level_data().get_expansion_length();
			int Lfrom = from.get_level_data().get_expansion_length();
			int L = std::max(Lto, Lfrom);
			cvector_t diag_coeffs(2 * L + 1, 1);
			for (int nu = -L; nu <= L; ++nu)
				diag_coeffs(nu + L) = std::exp(std::complex<double>(0., nu*(pi + theta)))
				* cyl_bessel_j(nu, z);
			return result_t(to.get_level_data(), from.get_level_data(),
				diag_coeffs);
		}
	};


	/** \brief the m2l operator */
	class m2l
		: public operator_with_wave_number<wave_number_t>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t cluster_t;
		typedef helmholtz_2d_wb_m2l_matrix result_t;

		m2l(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return m2l_indices<dimension>::eval(to.get_bounding_box(), from.get_bounding_box());
		}

		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			using boost::math::cyl_hankel_2;
			using namespace boost::math::double_constants;

			int L = to.get_level_data().get_expansion_length();
			location_t const &X = to.get_bounding_box().get_center();
			location_t const &Y = from.get_bounding_box().get_center();
			location_t d = X - Y;
			double r = d.norm();
			auto z = this->get_wave_number() * r;
			double theta = std::atan2(d(1), d(0));
			cvector_t diag_coeffs(2 * L + 1);
			for (int nu = -L; nu <= L; ++nu)
				diag_coeffs(nu + L) = std::exp(std::complex<double>(0., nu*(pi + theta)))
				* cyl_hankel_2(nu, z);
			return result_t(to.get_level_data(), diag_coeffs);
		}
	};


	/** \brief the p2m operator
	 * \tparam Ny the order of normal derivative
	 */
	template <unsigned int Ny>
	class p2m
		: public operator_with_wave_number<wave_number_t>
	{
	public:
		typedef operator_with_wave_number<wave_number_t> base_t;
		typedef helmholtz_2d_wb_fmm::cluster_t test_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_2d<>, WaveNumber>, 0, Ny
		>::trial_input_t trial_input_t;
		typedef cvector_t result_t;

		p2m(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		size_t rows(test_input_t const &to) const
		{
			return to.get_level_data().get_data_size();
		}

		result_t operator()(test_input_t const &to, trial_input_t const &y) const
		{
			result_t res = eval(to, y, std::integral_constant<int, Ny>());
			if (to.get_level_data().get_high_freq())
				return to.get_level_data().dft(res);
			return res;
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 0>) const
		{
			using boost::math::cyl_bessel_j;

			int L = to.get_level_data().get_expansion_length();
			location_t const &Y = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			location_t d = Y - y;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			result_t res(2 * L + 1, 1);
			auto z = this->get_wave_number() * r;
			for (int nu = -L; nu <= L; ++nu)
				res(-nu + L, 0) = std::exp(std::complex<double>(0, -nu * theta)) *
				cyl_bessel_j(nu, z);
			return res * std::complex<double>(0., -.25);
		}

		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 1>) const
		{
			using boost::math::cyl_bessel_j;
			using boost::math::cyl_bessel_j_prime;

			int L = to.get_level_data().get_expansion_length();
			location_t const &Y = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			location_t d = Y - y;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			double rdny = -d.dot(tri.get_unit_normal()) / r;
			location_t Td(d(1), -d(0));
			double thetadny = Td.dot(tri.get_unit_normal()) / (r * r);
			result_t res(2 * L + 1, 1);
			auto const &k = this->get_wave_number();
			auto z = k * r;
			for (int nu = -L; nu <= L; ++nu)
				res(-nu + L, 0) = std::exp(std::complex<double>(0, -nu * theta)) *
				(
					cyl_bessel_j_prime(nu, z) * k * rdny
					- std::complex<double>(0, nu * thetadny) * cyl_bessel_j(nu, z)
					);
			return res * std::complex<double>(0., -.25);
		}
	};

	/** \brief the p2l operator
	 * \tparam Ny the order of normal derivative
	 */
	template <unsigned int Ny>
	class p2l
		: public operator_with_wave_number<wave_number_t>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t test_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_2d<>, WaveNumber>, 0, Ny
		>::trial_input_t trial_input_t;
		typedef cvector_t result_t;

		p2l(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		size_t rows(test_input_t const &to) const
		{
			return to.get_level_data().get_data_size();
		}

		result_t operator()(test_input_t const &to, trial_input_t const &y) const
		{
			result_t res = eval(to, y, std::integral_constant<int, Ny>());
			if (to.get_level_data().get_high_freq())
				return to.get_level_data().dft(res);
			return res;
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 0>) const
		{
			using boost::math::cyl_hankel_2;
			using namespace boost::math::double_constants;

			int L = to.get_level_data().get_expansion_length();
			location_t const &X = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			location_t d = X - y;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			result_t res(2 * L + 1, 1);
			auto z = this->get_wave_number() * r;
			for (int nu = -L; nu <= L; ++nu)
				res(nu + L, 0) = std::exp(std::complex<double>(0., nu * (pi + theta))) *
				cyl_hankel_2(nu, z);
			return res * std::complex<double>(0., -.25);
		}

		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 1>) const
		{
			using boost::math::cyl_bessel_j_prime;
			using boost::math::cyl_neumann_prime;
			using boost::math::cyl_hankel_2;
			using namespace boost::math::double_constants;

			int L = to.get_level_data().get_expansion_length();
			location_t const &X = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			location_t d = X - y;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			double rdny = -d.dot(tri.get_unit_normal()) / r;
			location_t Td(d(1), -d(0));
			double thetadny = Td.dot(tri.get_unit_normal()) / (r * r);
			result_t res(2 * L + 1, 1);
			auto const &k = this->get_wave_number();
			auto z = k * r;
			for (int nu = -L; nu <= L; ++nu)
				res(nu + L, 0) =
				std::exp(std::complex<double>(0., nu * (pi + theta))) *
				(
				(cyl_bessel_j_prime(nu, z) - std::complex<double>(0, 1) * cyl_neumann_prime(nu, z)) * k * rdny +
					cyl_hankel_2(nu, z) * std::complex<double>(0, nu * thetadny)
					);
			return res * std::complex<double>(0., -.25);
		}
	};

	/** \brief the l2p operator
	 * \tparam Nx the order of normal derivative
	 */
	template <unsigned int Nx>
	class l2p
		: public operator_with_wave_number<wave_number_t>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t trial_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_2d<>, WaveNumber>, Nx, 0
		>::test_input_t test_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> result_t;

		l2p(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		size_t cols(trial_input_t const &from) const
		{
			return from.get_level_data().get_data_size();
		}

		result_t operator()(test_input_t const &x, trial_input_t const &from) const
		{
			result_t res = eval(x, from, std::integral_constant<int, Nx>());
			if (from.get_level_data().get_high_freq())
				return from.get_level_data().dft(res.conjugate().transpose()).conjugate().transpose() / double(cols(from));
			return res;
		}

	private:
		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 0>) const
		{
			using boost::math::cyl_bessel_j;

			int L = from.get_level_data().get_expansion_length();
			location_t const &X = from.get_bounding_box().get_center();
			location_t x = tsi.get_x();
			location_t d = x - X;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			auto z = this->get_wave_number() * r;
			result_t res(1, 2 * L + 1);
			for (int nu = -L; nu <= L; ++nu)
				res(0, nu + L) = std::exp(std::complex<double>(0, -nu * theta)) *
				cyl_bessel_j(nu, z);
			return res;
		}

		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 1>) const
		{
			using boost::math::cyl_bessel_j;
			using boost::math::cyl_bessel_j_prime;

			location_t const &X = from.get_bounding_box().get_center();
			location_t const &x = tsi.get_x();
			location_t d = x - X;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			double rdnx = d.dot(tsi.get_unit_normal()) / r;
			location_t Td(d(1), -d(0));
			double thetadnx = -Td.dot(tsi.get_unit_normal()) / (r * r);

			int L = from.get_level_data().get_expansion_length();
			result_t res(1, 2 * L + 1);

			auto const &k = this->get_wave_number();
			auto z = k * r;
			for (int nu = -L; nu <= L; ++nu)
				res(0, nu + L) = std::exp(std::complex<double>(0., -nu * theta)) *
				(
					cyl_bessel_j_prime(nu, z) * k * rdnx
					- cyl_bessel_j(nu, z) * std::complex<double>(0, nu * thetadnx)
					);
			return res;
		}
	};

	/** \brief the m2p operator
	 * \tparam Nx the order of normal derivative
	 */
	template <unsigned int Nx>
	class m2p
		: public operator_with_wave_number<wave_number_t>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t trial_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_2d<>, WaveNumber>, Nx, 0
		>::test_input_t test_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> result_t;

		m2p(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		size_t cols(trial_input_t const &from) const
		{
			return from.get_level_data().get_data_size();
		}

		result_t operator()(test_input_t const &x, trial_input_t const &from) const
		{
			result_t res = eval(x, from, std::integral_constant<int, Nx>());
			if (from.get_level_data().get_high_freq())
				return from.get_level_data().dft(res.conjugate().transpose()).conjugate().transpose() / double(cols(from));
			return res;
		}

	private:
		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 0>) const
		{
			using boost::math::cyl_hankel_2;
			using namespace boost::math::double_constants;

			int L = from.get_level_data().get_expansion_length();
			location_t const &Y = from.get_bounding_box().get_center();
			location_t const &x = tsi.get_x();
			auto d = x - Y;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			auto z = this->get_wave_number() * r;
			result_t res(1, 2 * L + 1);
			for (int nu = -L; nu <= L; ++nu)
				res(0, -nu + L) = std::exp(std::complex<double>(0, nu * (pi + theta))) *
				cyl_hankel_2(nu, z);
			return res;
		}

		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 1>) const
		{
			using boost::math::cyl_bessel_j_prime;
			using boost::math::cyl_neumann_prime;
			using boost::math::cyl_hankel_2;
			using namespace boost::math::double_constants;

			location_t const &Y = from.get_bounding_box().get_center();
			location_t const &x = tsi.get_x();
			location_t d = x - Y;
			double r = d.norm();
			double theta = std::atan2(d(1), d(0));
			double rdnx = d.dot(tsi.get_unit_normal()) / r;
			location_t Td(d(1), -d(0));
			double thetadnx = -Td.dot(tsi.get_unit_normal()) / (r * r);

			int L = from.get_level_data().get_expansion_length();
			result_t res(1, 2 * L + 1);

			auto const &k = this->get_wave_number();
			auto z = k * r;
			for (int nu = -L; nu <= L; ++nu)
				res(0, -nu + L) = std::exp(std::complex<double>(0., nu * (pi + theta))) *
				(
				(cyl_bessel_j_prime(nu, z) - std::complex<double>(0, 1) * cyl_neumann_prime(nu, z)) * k * rdnx
					+ cyl_hankel_2(nu, z) * std::complex<double>(0, nu * thetadnx)
					);
			return res;
		}
	};

	/** \brief constructor 
	 * \param [in] wave_number the wave number
	 */
	helmholtz_2d_wb_fmm(wave_number_t const &wave_number)
		: m_wave_number(wave_number)
	{
	}

	template <int Nx, int Ny>
	struct p2p_type
	{
		typedef fmm::p2p<NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_2d<>, wave_number_t>, Nx, Ny
		> > type;
	};

	template <int Ny>
	struct p2m_type
	{
		typedef p2m<Ny> type;
	};

	template <int Nx>
	struct m2p_type
	{
		typedef m2p<Nx> type;
	};

	template <int Ny>
	struct p2l_type
	{
		typedef p2l<Ny> type;
	};

	template <int Nx>
	struct l2p_type
	{
		typedef l2p<Nx> type;
	};

	template <int Nx, int Ny>
	fmm::p2p<NiHu::normal_derivative_kernel<
		NiHu::helmholtz_kernel<NiHu::space_2d<>, wave_number_t>, Nx, Ny
	> >
		create_p2p() const
	{
		typedef NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_2d<>, wave_number_t>, Nx, Ny
		> kernel_t;
		return fmm::p2p<kernel_t>(kernel_t(m_wave_number));
	}

	template <int Ny>
	p2m<Ny> create_p2m() const
	{
		return p2m<Ny>(this->get_wave_number());
	}

	template <int Ny>
	p2l<Ny> create_p2l() const
	{
		return p2l<Ny>(this->get_wave_number());
	}

	template <int Nx>
	m2p<Nx> create_m2p() const
	{
		return m2p<Nx>(this->get_wave_number());
	}

	template <int Nx>
	l2p<Nx> create_l2p() const
	{
		return l2p<Nx>(this->get_wave_number());
	}

	m2m create_m2m() const
	{
		return m2m(this->get_wave_number());
	}

	l2l create_l2l() const
	{
		return l2l(this->get_wave_number());
	}

	m2l create_m2l() const
	{
		return m2l(this->get_wave_number());
	}

	wave_number_t const &get_wave_number() const
	{
		return m_wave_number;
	}

	void set_accuracy(double)
	{

	}

	void init_level_data(cluster_tree_t const &tree)
	{
		/// \todo make sure that leaf level is in low frequency domain

		using namespace boost::math::double_constants;
		double lambda = two_pi / std::real(m_wave_number);
		m_level_data_vector.clear();
		m_level_data_vector.resize(tree.get_n_levels());
		for (size_t i = 2; i < tree.get_n_levels(); ++i)
		{
			size_t idx = tree.level_begin(i);
			double D = tree[idx].get_bounding_box().get_diameter();
			m_level_data_vector[i].init(D / lambda);
		}

		for (size_t i = 2; i < tree.get_n_levels() - 1; ++i)
		{
			m_level_data_vector[i].set_interp_up(
				m_level_data_vector[i + 1].get_expansion_length());
		}

		for (size_t i = 3; i < tree.get_n_levels(); ++i)
		{
			m_level_data_vector[i].set_interp_dn(
				m_level_data_vector[i - 1].get_expansion_length());
		}
	}

	helmholtz_2d_wb_level_data const &get_level_data(size_t idx) const
	{
		return m_level_data_vector[idx];
	}

	helmholtz_2d_wb_level_data &get_level_data(size_t idx)
	{
		return m_level_data_vector[idx];
	}

	void print_level_data(std::ostream &os = std::cout) const
	{
		for (size_t i = 2; i < m_level_data_vector.size(); ++i)
		{
			auto const &ld = m_level_data_vector[i];
			os << "level: " << i << " " << ld.get_expansion_length() << ' ' << ld.get_high_freq() << std::endl;
		}
	}

private:
	wave_number_t m_wave_number;
	std::vector<helmholtz_2d_wb_level_data> m_level_data_vector;
};

} // end of namespace fmm
} // namespace NiHu

#endif
