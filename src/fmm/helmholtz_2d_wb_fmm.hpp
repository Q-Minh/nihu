/** 
 * \file helmholtz_2d_wb_fmm.hpp
 * \brief 2D Wide Band Helmholtz FMM
 * \ingroup fmm_helmholtz_2d_wb
 */

#ifndef HELMHOLTZ_2D_WB_FMM_HPP_INCLUDED
#define HELMHOLTZ_2D_WB_FMM_HPP_INCLUDED

#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
#include <boost/math/constants/constants.hpp>

#include "cluster.hpp"
#include "cluster_tree.hpp"
#include "fmm_operator.hpp"
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
	/// \brief the cluster type of the FMM
	typedef helmholtz_2d_wb_cluster cluster_t;
	/// \brief the cluster tree type
	typedef cluster_tree<cluster_t> cluster_tree_t;
	/// \brief a complex vector type
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;
	/// \brief the space dimension
	static size_t const dimension = cluster_t::dimension;
	/// \brief the cluster's bounding box type
	typedef typename cluster_t::bounding_box_t bounding_box_t;
	/// \brief the physical location type
	typedef typename bounding_box_t::location_t location_t;

private:
	typedef NiHu::helmholtz_kernel<NiHu::space_2d<>, wave_number_t> distance_dependent_kernel_t;

private:
	/// \brief convert from 2D Cartesian to polar coordinates
	/// \param [in] d the location vector
	/// \param [out] r the distance
	/// \param [out] theta the angle
	static void cart2pol(location_t const &d, double &r, double &theta)
	{
		r = d.norm();
		theta = std::atan2(d(1), d(0));
	}

public:
	/** \brief the m2m operator */
	class m2m
		: public operator_with_wave_number<wave_number_t>
		, public fmm_operator<m2m_tag>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		/// \brief the cluster type
		typedef helmholtz_2d_wb_fmm::cluster_t cluster_t;
		/// \brief the evaluated operator type
		typedef helmholtz_2d_wb_m2m_matrix result_t;

		/// \brief constructor of the operator
		/// \param [in] wave_number the wave number
		m2m(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		/// \brief return a unique index for a source and receiver cluster
		/// \param [in] to the receiver cluster
		/// \param [in] from the source cluster
		/// \return a unique operator index
		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return bounding_box_t::dist2idx(
				from.get_bounding_box().get_center(),
				to.get_bounding_box().get_center());
		}

	public:

		/// \brief evaluate the operator for a source and receiver cluster
		/// \param [in] to the receiver cluster
		/// \param [in] from the source cluster
		/// \return the evaluated operator
		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			using boost::math::cyl_bessel_j;

			location_t const &X = to.get_bounding_box().get_center();
			location_t const &Y = from.get_bounding_box().get_center();
			double r, theta;
			cart2pol(X - Y, r, theta);
			auto z = this->get_wave_number() * r;
			int Lto = int(to.get_level_data().get_expansion_length());
			int Lfrom = int(from.get_level_data().get_expansion_length());
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
		, public fmm_operator<l2l_tag>
	{
	private:
		typedef operator_with_wave_number<WaveNumber> base_t;

	public:
		/// \brief the cluster type
		typedef helmholtz_2d_wb_fmm::cluster_t cluster_t;
		/// \brief the evaluated operator's type
		typedef helmholtz_2d_wb_l2l_matrix result_t;

		/// \brief constructor
		/// \param [in] wave_number the wave number
		l2l(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		/// \brief return a unique index for a source and receiver cluster
		/// \param [in] to receiver cluster
		/// \param [in] from source cluster
		/// \return the unique operator index for this cluster pair
		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return bounding_box_t::dist2idx(
				to.get_bounding_box().get_center(),
				from.get_bounding_box().get_center());
		}

		/// \brief evaluate the operator for a source and receiver cluster
		/// \param [in] to the receiver cluster
		/// \param [in] from the source cluster
		/// \return the evaluated operator
		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			using boost::math::cyl_bessel_j;
			using namespace boost::math::double_constants;

			location_t const &X = to.get_bounding_box().get_center();
			location_t const &Y = from.get_bounding_box().get_center();
			double r, theta;
			cart2pol(X - Y, r, theta);
			auto z = this->get_wave_number() * r;
			int Lto = int(to.get_level_data().get_expansion_length());
			int Lfrom = int(from.get_level_data().get_expansion_length());
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
		, public fmm_operator<m2l_tag>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		/// \brief the cluster type
		typedef helmholtz_2d_wb_fmm::cluster_t cluster_t;
		/// \the evaluated operator's type
		typedef helmholtz_2d_wb_m2l_matrix result_t;

		/// \brief constructor
		/// \param [in] wave_number the wave number
		m2l(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		/// \brief assign a unique index to source and receiver clusters
		/// \param [in] to the receiver cluster
		/// \param [in] from the source cluster
		/// \return a unique operator index
		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return m2l_indices<dimension>::eval(to.get_bounding_box(), from.get_bounding_box());
		}

		/// \brief evaluate the operator for a source and receiver cluster
		/// \param [in] to the receiver cluster
		/// \param [in] from the source cluster
		/// \return the evaluated operator
		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			using boost::math::cyl_hankel_2;
			using namespace boost::math::double_constants;

			size_t L = to.get_level_data().get_expansion_length();
			location_t const &X = to.get_bounding_box().get_center();
			location_t const &Y = from.get_bounding_box().get_center();
			double r, theta;
			cart2pol(X - Y, r, theta);
			auto z = this->get_wave_number() * r;
			cvector_t diag_coeffs(2 * L + 1);
			for (int nu = -int(L); nu <= int(L); ++nu)
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
		, public fmm_operator<p2m_tag>
	{
	public:
		typedef operator_with_wave_number<wave_number_t> base_t;
		/// \brief the test input type
		typedef helmholtz_2d_wb_fmm::cluster_t test_input_t;
		/// \brief the trial input type
		typedef typename NiHu::normal_derivative_kernel<
			distance_dependent_kernel_t, 0, Ny
		>::trial_input_t trial_input_t;
		/// \brief the evaluated operator's type
		typedef cvector_t result_t;

		/// \brief constructor of the operator
		/// \param [in] wave_number the wave number
		p2m(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		/// \brief number of rows of the operator
		/// \param [in] to the receiver
		/// \return number of rows
		size_t rows(test_input_t const &to) const
		{
			return to.get_level_data().get_data_size();
		}

		/// \brief evaluate the operator for a source and receiver
		/// \param [in] to the receiver
		/// \param [in] from the source
		/// \return the evaluated operator
		result_t operator()(test_input_t const &to, trial_input_t const &y) const
		{
			result_t res = eval(to, y, std::integral_constant<int, Ny>());
			if (to.get_level_data().is_high_freq())
				return to.get_level_data().dft(res);
			return res;
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 0>) const
		{
			using boost::math::cyl_bessel_j;

			int L = int(to.get_level_data().get_expansion_length());
			location_t const &Y = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			double r, theta;
			cart2pol(Y - y, r, theta);
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

			int L = int(to.get_level_data().get_expansion_length());
			location_t const &Y = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			location_t d = Y - y;
			double r, theta;
			cart2pol(d, r, theta);
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
		, public fmm_operator<p2l_tag>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t test_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			distance_dependent_kernel_t, 0, Ny
		>::trial_input_t trial_input_t;
		typedef cvector_t result_t;

		/// \brief constructor of the operator
		/// \param [in] wave_number the wave number
		p2l(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		/// \brief number of rows of the operator
		/// \param [in] to the receiver
		/// \return number of rows
		size_t rows(test_input_t const &to) const
		{
			return to.get_level_data().get_data_size();
		}

		/// \brief evaluate the operator for a source and receiver
		/// \param [in] to the receiver
		/// \param [in] y the source
		/// \return the evaluated operator
		result_t operator()(test_input_t const &to, trial_input_t const &y) const
		{
			result_t res = eval(to, y, std::integral_constant<int, Ny>());
			if (to.get_level_data().is_high_freq())
				return to.get_level_data().dft(res);
			return res;
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 0>) const
		{
			using boost::math::cyl_hankel_2;
			using namespace boost::math::double_constants;

			int L = int(to.get_level_data().get_expansion_length());
			location_t const &X = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			double r, theta;
			cart2pol(X - y, r, theta);
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

			int L = int(to.get_level_data().get_expansion_length());
			location_t const &X = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			location_t d = X - y;
			double r, theta;
			cart2pol(X - y, r, theta);
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
		, public fmm_operator<l2p_tag>

	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t trial_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			distance_dependent_kernel_t, Nx, 0
		>::test_input_t test_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> result_t;

		/// \brief constructor of the operator
		/// \param [in] wave_number the wave number
		l2p(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		/// \brief number of columns of the operator
		/// \param [in] from the source
		/// \return number of columns
		size_t cols(trial_input_t const &from) const
		{
			return from.get_level_data().get_data_size();
		}

		/// \brief evaluate the operator for a source and receiver
		/// \param [in] x the receiver
		/// \param [in] from the source
		/// \return the evaluated operator
		result_t operator()(test_input_t const &x, trial_input_t const &from) const
		{
			result_t res = eval(x, from, std::integral_constant<int, Nx>());
			if (from.get_level_data().is_high_freq())
				return from.get_level_data().dft(res.conjugate().transpose()).conjugate().transpose() / double(cols(from));
			return res;
		}

	private:
		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 0>) const
		{
			using boost::math::cyl_bessel_j;

			int L = int(from.get_level_data().get_expansion_length());
			location_t const &X = from.get_bounding_box().get_center();
			location_t x = tsi.get_x();
			double r, theta;
			cart2pol(x - X, r, theta);
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
			double r, theta;
			cart2pol(d, r, theta);
			double rdnx = d.dot(tsi.get_unit_normal()) / r;
			location_t Td(d(1), -d(0));
			double thetadnx = -Td.dot(tsi.get_unit_normal()) / (r * r);

			int L = int(from.get_level_data().get_expansion_length());
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
		, public fmm_operator<m2p_tag>
	{
	private:
		typedef operator_with_wave_number<wave_number_t> base_t;

	public:
		typedef helmholtz_2d_wb_fmm::cluster_t trial_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			distance_dependent_kernel_t, Nx, 0
		>::test_input_t test_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> result_t;

		/// \brief constructor of the operator
		/// \param [in] wave_number the wave number
		m2p(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		/// \brief number of columns of the operator
		/// \param [in] from the source
		/// \return number of columns
		size_t cols(trial_input_t const &from) const
		{
			return from.get_level_data().get_data_size();
		}

		/// \brief evaluate the operator for a source and receiver
		/// \param [in] x the receiver
		/// \param [in] from the source
		/// \return the evaluated operator
		result_t operator()(test_input_t const &x, trial_input_t const &from) const
		{
			result_t res = eval(x, from, std::integral_constant<int, Nx>());
			if (from.get_level_data().is_high_freq())
				return from.get_level_data().dft(res.conjugate().transpose()).conjugate().transpose() / double(cols(from));
			return res;
		}

	private:
		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 0>) const
		{
			using boost::math::cyl_hankel_2;
			using namespace boost::math::double_constants;

			int L = int(from.get_level_data().get_expansion_length());
			location_t const &Y = from.get_bounding_box().get_center();
			location_t const &x = tsi.get_x();
			auto d = x - Y;
			double r, theta;
			cart2pol(d, r, theta);
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
			double r, theta;
			cart2pol(d, r, theta);
			double rdnx = d.dot(tsi.get_unit_normal()) / r;
			location_t Td(d(1), -d(0));
			double thetadnx = -Td.dot(tsi.get_unit_normal()) / (r * r);

			int L = int(from.get_level_data().get_expansion_length());
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
			distance_dependent_kernel_t, Nx, Ny
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
		distance_dependent_kernel_t, Nx, Ny
	> >
		create_p2p() const
	{
		typedef NiHu::normal_derivative_kernel<
			distance_dependent_kernel_t, Nx, Ny> kernel_t;
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

	/// \brief set the method's prescribed accuracy
	void set_accuracy(double)
	{
	}

	/// \brief initialize level data for a specified cluster tree
	/// \param [in] tree the cluster tree
	void init_level_data(cluster_tree_t const &tree)
	{
		using namespace boost::math::double_constants;

		double lambda = two_pi / std::real(m_wave_number);
		m_level_data_vector.clear();
		m_level_data_vector.resize(tree.get_n_levels());

		// initialize level data for each level where anything happens
		for (size_t i = 2; i < tree.get_n_levels(); ++i)
		{
			// index of first cluster at the level
			size_t idx = tree.level_begin(i);
			double D = tree[idx].get_bounding_box().get_diameter();
			m_level_data_vector[i].init(D / lambda);
		}

		// check that the highest level leaf cluster is in the low frequency domain
		size_t idx;
		for (idx = 0; idx < tree.get_n_clusters(); ++idx)
			if (tree[idx].is_leaf() && m_level_data_vector[tree[idx].get_level()].is_high_freq())
				throw std::runtime_error("Leaf level cluster (" + std::to_string(idx) + ") has been found in the high frequency domain.");

		// initialize upward interpolation for each m2m receiver level
		for (size_t i = 2; i < tree.get_n_levels() - 1; ++i)
		{
			size_t L_from = m_level_data_vector[i + 1].get_expansion_length();
			m_level_data_vector[i].set_interp_up(L_from);
		}

		// initialize downward interpolation for each l2l receiver level
		for (size_t i = 3; i < tree.get_n_levels(); ++i)
		{
			size_t L_from = m_level_data_vector[i - 1].get_expansion_length();
			m_level_data_vector[i].set_interp_dn(L_from);
		}
	}

	/// \brief return level data at a specified level
	/// \param [in] idx the level index
	/// \return level data
	helmholtz_2d_wb_level_data const &get_level_data(size_t idx) const
	{
		return m_level_data_vector[idx];
	}

	/// \brief return level data reference at a specified level
	/// \param [in] idx the level index
	/// \return level data reference
	helmholtz_2d_wb_level_data &get_level_data(size_t idx)
	{
		return m_level_data_vector[idx];
	}

	/// \brief print debug information
	/// \param [in] os the output stream
	void print_level_data(std::ostream &os = std::cout) const
	{
		for (size_t i = 2; i < m_level_data_vector.size(); ++i)
		{
			auto const &ld = m_level_data_vector[i];
			os << "level: " << i << " " << ld.get_expansion_length() << ' ' << ld.is_high_freq() << std::endl;
		}
	}

private:
	wave_number_t m_wave_number;
	std::vector<helmholtz_2d_wb_level_data> m_level_data_vector;
};

} // end of namespace fmm
} // namespace NiHu

#endif
