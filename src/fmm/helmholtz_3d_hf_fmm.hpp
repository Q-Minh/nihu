/** \file helmholtz_3d_hf_fmm.hpp
 * \brief implemenetation of the 3D High Frequency Helmholtz FMM
 */
#ifndef HELMHOLTZ_3D_HF_FMM_HPP_INCLUDED
#define HELMHOLTZ_3D_HF_FMM_HPP_INCLUDED

#include "cluster.hpp"
#include "cluster_tree.hpp"
#include "helmholtz_3d_hf_cluster.h"
#include "helmholtz_3d_hf_level_data.h"
#include "m2l_indices.hpp"
#include "operator_with_wave_number.hpp"
#include "p2p.hpp"
#include "unit_sphere.h"

#include "library/helmholtz_kernel.hpp"
#include "library/helmholtz_singular_integrals.hpp"
#include "library/helmholtz_nearly_singular_integrals.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <Eigen/Dense>
#include <complex>

namespace NiHu
{
namespace fmm
{


class up_shift
{
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

public:
	up_shift()
		: m_level_data(nullptr)
	{
	}

	up_shift(cvector_t const &shift, helmholtz_3d_hf_level_data const &ld)
		: m_shift(shift)
		, m_level_data(&ld)
	{
	}

	/// \brief multiply the updownshift matrix with a multipole / local
	/// \param rhs the right hand side
	cvector_t operator*(cvector_t const &rhs) const
	{
		return m_shift.array() * m_level_data->interp_up(rhs).array();
	}

private:
	cvector_t m_shift;
	helmholtz_3d_hf_level_data const *m_level_data;
};


class down_shift
{
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

public:
	down_shift()
		: m_level_data(nullptr)
	{
	}

	down_shift(cvector_t const &shift, helmholtz_3d_hf_level_data const &ld)
		: m_shift(shift)
		, m_level_data(&ld)
	{
	}

	/// \brief multiply the updownshift matrix with a multipole / local
	/// \param [in] rhs the right hand side
	cvector_t operator*(cvector_t const &rhs) const
	{
		return m_shift.array() * m_level_data->interp_down(rhs).array();
	}

private:
	cvector_t m_shift;
	helmholtz_3d_hf_level_data const *m_level_data;
};


/// \brief the fmm for the 3D Helmholtz equation
/// \tparam WaveNumber the wave number type (double or complex)
template <class WaveNumber>
class helmholtz_3d_hf_fmm
{
public:
	typedef WaveNumber wave_number_t;

	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

	typedef helmholtz_3d_hf_cluster cluster_t;
	typedef typename cluster_t::bounding_box_t bounding_box_t;
	typedef typename bounding_box_t::location_t location_t;
	typedef cluster_tree<cluster_t> cluster_tree_t;

	helmholtz_3d_hf_fmm(wave_number_t const &k)
		: m_wave_number(k)
	{
	}

private:
	static size_t compute_expansion_length(double drel, double C)
	{
		double const pi = boost::math::constants::pi<double>();
		double kd = 2. * pi * drel;
		return size_t(std::ceil(kd + C * std::log(kd + pi)));
	}

public:
	void set_accuracy(double C)
	{
		m_C = C;
	}

	void init_level_data(cluster_tree_t const &tree)
	{
		double const pi = boost::math::constants::pi<double>();
		double lambda = 2. * pi / std::real(m_wave_number);
		m_level_data_vector.clear();
		m_level_data_vector.resize(tree.get_n_levels());
		for (size_t i = 0; i < tree.get_n_levels(); ++i)
		{
			auto &ld = m_level_data_vector[i];
			size_t idx = tree.level_begin(i);
			double d = tree[idx].get_bounding_box().get_diameter();
			ld.set_expansion_length(compute_expansion_length(d / lambda, m_C));
		}

		for (size_t i = 0; i < tree.get_n_levels(); ++i)
		{
			auto &ld = m_level_data_vector[i];
			auto const &Sto = ld.get_unit_sphere();
			std::cout << "--- now --- " << std::endl;
			if (i != 0)
				ld.set_interp_dn(interpolator(m_level_data_vector[i - 1].get_unit_sphere(), Sto));
			if (i != tree.get_n_levels() - 1)
				ld.set_interp_up(interpolator(m_level_data_vector[i + 1].get_unit_sphere(), Sto));
		}
	}

	helmholtz_3d_hf_level_data const &get_level_data(size_t idx) const
	{
		return m_level_data_vector[idx];
	}

	helmholtz_3d_hf_level_data &get_level_data(size_t idx)
	{
		return m_level_data_vector[idx];
	}

	/// \brief M2M operator of the FMM for the Helmholtz equation in 3D
	class m2m
		: public operator_with_wave_number<wave_number_t>
	{
	public:
		typedef operator_with_wave_number<wave_number_t> base_t;
		typedef helmholtz_3d_hf_fmm::cluster_t cluster_t;
		typedef up_shift result_t;

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
			auto const &Sto = to.get_level_data().get_unit_sphere();
			location_t D = to.get_bounding_box().get_center()
				- from.get_bounding_box().get_center();
			auto const &k = this->get_wave_number();
			std::complex<double> const J(0., 1.);
			cvector_t shift = exp((-J * k * Sto.get_s().transpose() * D).array());
			return result_t(shift, to.get_level_data());
		}
	};


	/// \brief L2L operator of the FMM for the Helmholtz equation in 3D
	class l2l
		: public operator_with_wave_number<wave_number_t>
	{
	public:
		typedef operator_with_wave_number<WaveNumber> base_t;
		typedef helmholtz_3d_hf_fmm::cluster_t cluster_t;
		typedef down_shift result_t;

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
			auto const &Sto = to.get_level_data().get_unit_sphere();
			location_t D = to.get_bounding_box().get_center()
				- from.get_bounding_box().get_center();
			auto const &k = this->get_wave_number();
			std::complex<double> const J(0., 1.);
			cvector_t shift = exp((-J * k * Sto.get_s().transpose() * D).array());
			return result_t(shift, to.get_level_data());
		}
	};


	/// \brief P2M operator of the FMM for the Helmholtz equation in 3D
	/// \tparam Ny the order of differentiation w.r.t y
	template <unsigned int Ny>
	class p2m
		: public operator_with_wave_number<wave_number_t>
	{
	public:
		typedef operator_with_wave_number<wave_number_t> base_t;
		typedef typename cluster_t::location_t location_t;

		typedef cluster_t test_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_3d<>, WaveNumber>, 0, Ny
		>::trial_input_t trial_input_t;
		typedef cvector_t result_t;

		p2m(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		size_t rows(test_input_t const &to) const
		{
			return to.get_level_data().get_unit_sphere().get_s().cols();
		}

		result_t operator()(test_input_t const &to, trial_input_t const &tri) const
		{
			return eval(to, tri, std::integral_constant<int, Ny>());
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 0>) const
		{
			location_t const &Y = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			auto const &s = to.get_level_data().get_unit_sphere().get_s();
			auto const &k = this->get_wave_number();
			std::complex<double> const J(0.0, 1.0);
			location_t d = Y - y;
			return Eigen::exp(-J * k*(s.transpose() * d).array());
		}

		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 1>) const
		{
			location_t const &Y = to.get_bounding_box().get_center();
			location_t const &y = tri.get_x();
			auto const &s = to.get_level_data().get_unit_sphere().get_s();
			auto const &k = this->get_wave_number();
			std::complex<double> const J(0.0, 1.0);
			location_t d = Y - y;
			return Eigen::exp(-J * k*(s.transpose() * d).array())
				* ((J*k) * (s.transpose() * tri.get_unit_normal()).array());
		}
	};

	/// \brief P2L operator of the FMM for the Helmholtz equation in 3D
	/// \tparam Ny the order of differentiation w.r.t y
	template <unsigned int Ny>
	class p2l
		: public operator_with_wave_number<wave_number_t>
	{
	public:
		typedef operator_with_wave_number<wave_number_t> base_t;

		typedef typename cluster_t::location_t location_t;

		typedef cluster_t test_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_3d<>, WaveNumber>, 0, Ny
		>::trial_input_t trial_input_t;
		typedef cvector_t result_t;

		p2l(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		size_t rows(test_input_t const &to) const
		{
			return to.get_level_data().get_unit_sphere().get_s().cols();
		}

		result_t operator()(test_input_t const &to, trial_input_t const &tri) const
		{
			return eval(to, tri, std::integral_constant<int, Ny>());
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 0>) const
		{
			throw std::logic_error("Unimplemented p2l operator");
			return result_t::Zero(rows(to), 1);
		}

		result_t eval(test_input_t const &to, trial_input_t const &tri,
			std::integral_constant<int, 1>) const
		{
			throw std::logic_error("Unimplemented p2l operator");
			return result_t::Zero(rows(to), 1);
		}
	};


	/// \brief L2P operator of the FMM for the Helmholtz equation in 3D
	/// \tparam Nx the order of differentiation w.r.t x
	template <unsigned int Nx>
	class l2p
		: public operator_with_wave_number<wave_number_t>
	{
	public:
		typedef operator_with_wave_number<wave_number_t> base_t;

		typedef typename cluster_t::location_t location_t;

		typedef cluster_t trial_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_3d<>, WaveNumber>, Nx, 0
		>::test_input_t test_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> result_t;

		l2p(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		size_t cols(trial_input_t const &from) const
		{
			return from.get_level_data().get_unit_sphere().get_s().cols();
		}

		result_t operator()(test_input_t const &tsi, trial_input_t const &from) const
		{
			return eval(tsi, from, std::integral_constant<int, Nx>());
		}

	private:
		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 0>) const
		{
			location_t const &X = from.get_bounding_box().get_center();
			location_t const &x = tsi.get_x();
			auto const &S = from.get_level_data().get_unit_sphere();
			auto const &s = S.get_s();
			auto const &w = S.get_w();
			auto const &k = this->get_wave_number();
			std::complex<double> const J(0.0, 1.0);
			location_t d = x - X;
			return Eigen::exp(-J * k*(s.transpose() * d).array()) * w.array();
		}

		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 1>) const
		{
			location_t const &X = from.get_bounding_box().get_center();
			location_t const &x = tsi.get_x();
			auto const &S = from.get_level_data().get_unit_sphere();
			auto const &s = S.get_s();
			auto const &w = S.get_w();
			auto const &k = this->get_wave_number();
			std::complex<double> const J(0.0, 1.0);
			location_t d = x - X;
			return Eigen::exp(-J * k*(s.transpose() * d).array())
				* (-J * k) * (s.transpose() * tsi.get_unit_normal()).array()
				* w.array();
		}
	};


	/// \brief M2P operator of the FMM for the Helmholtz equation in 3D
	/// \tparam Nx the order of differentiation w.r.t x
	template <unsigned int Nx>
	class m2p
		: public operator_with_wave_number<wave_number_t>
	{
	public:
		typedef operator_with_wave_number<wave_number_t> base_t;

		typedef typename cluster_t::location_t location_t;

		typedef cluster_t trial_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_3d<>, WaveNumber>, Nx, 0
		>::test_input_t test_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> result_t;

		m2p(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		size_t cols(trial_input_t const &from) const
		{
			return from.get_level_data().get_unit_sphere().get_s().cols();
		}

		result_t operator()(test_input_t const &tsi, trial_input_t const &from) const
		{
			return eval(tsi, from, std::integral_constant<int, Nx>());
		}

	private:
		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 0>) const
		{
			throw std::logic_error("Unimplemented m2p operator");
			return result_t::Zero(1, cols(from));
		}

		result_t eval(test_input_t const &tsi, trial_input_t const &from,
			std::integral_constant<int, 1>) const
		{
			throw std::logic_error("Unimplemented m2p operator");
			return result_t::Zero(1, cols(from));
		}
	};

	/// \brief M2L operator of the FMM for the Helmholtz equation in 3D
	class m2l
		: public operator_with_wave_number<wave_number_t>
	{
	public:
		typedef operator_with_wave_number<wave_number_t> base_t;
		typedef helmholtz_3d_hf_fmm::cluster_t cluster_t;
		typedef Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic> result_t;

		m2l(wave_number_t const &wave_number)
			: base_t(wave_number)
		{
		}

		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return m2l_indices<3>::eval(to.get_bounding_box(), from.get_bounding_box());
		}

		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			auto const &k = this->get_wave_number();

			location_t const &X = to.get_bounding_box().get_center();
			location_t const &Y = from.get_bounding_box().get_center();
			location_t Dvec = X - Y;

			size_t L = to.get_level_data().get_expansion_length();
			auto const &s = to.get_level_data().get_unit_sphere().get_s();

			return m2l_matrix_impl(Dvec, s, k, L).asDiagonal();
		}

	private:
		static cvector_t m2l_matrix_impl(location_t const &Dvec,
				Eigen::Matrix<double, 3, Eigen::Dynamic> const &s,
				wave_number_t const &k,
				size_t L)
		{
			double const pi = boost::math::constants::pi<double>();
			double D = Dvec.norm();
			location_t Uvec = Dvec / D;
			size_t N = s.cols();

			Eigen::Matrix<double, 1, Eigen::Dynamic> x = Uvec.transpose() * s;
			auto z = -k * D;

			cvector_t M2L = cvector_t::Zero(N, 1);
			std::complex<double> const J(0.0, 1.0);

			for (size_t l = 0; l <= L; ++l)
			{
				auto h = boost::math::sph_hankel_1(l, z);
				auto c = double(2 * l + 1) * std::pow(J, l) * h;
				for (int i = 0; i < s.cols(); ++i)
					M2L(i) += c * boost::math::legendre_p(int(l), x(0, i));
			}
			M2L *= -J * k / (4.0*pi) / (4.0*pi);

			return M2L;
		}
	};


	template <int Ny>
	p2m<Ny> create_p2m() const
	{
		return p2m<Ny>(m_wave_number);
	}

	template <int Ny>
	p2l<Ny> create_p2l() const
	{
		return p2l<Ny>(m_wave_number);
	}

	template <int Nx>
	l2p<Nx> create_l2p() const
	{
		return l2p<Nx>(m_wave_number);
	}

	template <int Nx>
	m2p<Nx> create_m2p() const
	{
		return m2p<Nx>(m_wave_number);
	}

	template <int Nx, int Ny>
	fmm::p2p<NiHu::normal_derivative_kernel<
		NiHu::helmholtz_kernel<NiHu::space_3d<>, WaveNumber>, Nx, Ny
		> > create_p2p() const
	{
		typedef NiHu::normal_derivative_kernel<
			NiHu::helmholtz_kernel<NiHu::space_3d<>, WaveNumber>, Nx, Ny
		> kernel_t;
		return fmm::p2p<kernel_t>(kernel_t(m_wave_number));
	}

	m2m create_m2m() const
	{
		return m2m(m_wave_number);
	}

	l2l create_l2l() const
	{
		return l2l(m_wave_number);
	}

	m2l create_m2l() const
	{
		return m2l(m_wave_number);
	}

private:
	wave_number_t m_wave_number;
	std::vector<helmholtz_3d_hf_level_data> m_level_data_vector;
	double m_C;
};


} // end of namespace fmm
} // namespace NiHu

#endif // HELMHOLTZ_3D_HF_FMM_HPP_INCLUDED
