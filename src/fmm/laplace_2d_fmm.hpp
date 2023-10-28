/** \file laplace_2d_fmm.hpp
 * \brief implementation of the FMM method for the 2D Laplace equation
 */
#ifndef LAPLACE_2D_FMM_HPP_INCLUDED
#define LAPLACE_2D_FMM_HPP_INCLUDED

#include <boost/math/constants/constants.hpp>

#include "cluster.hpp"
#include "m2l_indices.hpp"
#include "p2p.hpp"

#include "library/laplace_kernel.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/laplace_singular_integrals.hpp"

#include <boost/math/special_functions/binomial.hpp>

#include <complex>

namespace NiHu
{
namespace fmm
{

/// \brief cluster type of the Laplace 2D FMM
class laplace_2d_cluster;

/// \brief specialisation of cluster traits for the 2D Laplace FMM
template <>
struct cluster_traits<laplace_2d_cluster>
{
	/// \brief the space dimension
	static size_t const dimension = 2;
	/// \brief the multipole type
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> multipole_t;
	/// \brief the local type
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> local_t;
};

/// \brief cluster type of the Laplace 2D FMM
class laplace_2d_cluster
	: public cluster_base<laplace_2d_cluster>
{
	typedef cluster_base<laplace_2d_cluster> base_t;
	
public:
	/// \brief the multipole type
	typedef typename base_t::multipole_t multipole_t;
	/// \brief the local type
	typedef typename base_t::local_t local_t;

	/// \brief set the expansion length
	void set_expansion_length(size_t L);

	/// \brief get the expansion length
	/// \return the expansion length
	size_t get_expansion_length() const;

	/// \brief return a cleared multipole contribution
	multipole_t zero_multipole() const;

	/// \brief return a cleared local contribution
	local_t zero_local() const;

private:
	size_t m_expansion_length;
};

/// \brief the fmm for the Laplace equation in 2D
class laplace_2d_fmm
{
public:
	/// \brief the cluster type
	typedef laplace_2d_cluster cluster_t;
	/// \brief the bounding box type
	typedef cluster_t::bounding_box_t bounding_box_t;
	/// \brief the location type
	typedef cluster_t::location_t location_t;

private:
	/// \brief helper function to convert a cluster's center into a complex number
	/// \param [in] c the cluster
	/// \return the cluster's center interpreted as a complex number
	static std::complex<double> center2complex(cluster_t const &c);

public:
	/// \brief return an instance of the P2P operator
	/// \tparam Nx order of differentiation w.r.t x
	/// \tparam Ny order of differentiation w.r.t y
	/// \return an instance of the P2P operator
	template <int Nx, int Ny>
	fmm::p2p<NiHu::normal_derivative_kernel<
		NiHu::laplace_kernel<NiHu::space_2d<double> >, Nx, Ny
	> >
		create_p2p() const
	{
		typedef NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_2d<double> >, Nx, Ny
		> kernel_t;
		return fmm::p2p<kernel_t>(kernel_t());
	}


	/// \brief P2M operator of the Laplace 2D FMM
	/// \tparam Ny order of derivation
	template <unsigned Ny>
	class p2m
	{
	public:
		typedef cluster_t test_input_t;

		typedef typename NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_2d<double> >,
			0, Ny
		>::trial_input_t trial_input_t;

		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> result_t;

		size_t rows(test_input_t const &to) const
		{
			return to.get_expansion_length() + 1;
		}

		result_t operator()(test_input_t const &to, trial_input_t const & from) const
		{
			return eval(to, from, std::integral_constant<unsigned, Ny>());
		}

		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 0>) const
		{
			using namespace boost::math::double_constants;
			
			auto const& y = from.get_x();
			std::complex<double> zy(y(0), y(1));
			zy -= center2complex(to);

			result_t res(to.get_expansion_length() + 1, 1);
			res(0) = 1.0;
			for (size_t k = 1; k <= to.get_expansion_length(); ++k)
				res(k) = -std::pow(zy, k) / double(k);
			return -res / two_pi;
		}

		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 1>) const
		{
			using namespace boost::math::double_constants;
			
			auto const& y = from.get_x();
			std::complex<double> zy(y(0), y(1));
			zy -= center2complex(to);
			std::complex<double> zydny(from.get_unit_normal()(0), from.get_unit_normal()(1));
			result_t res(to.get_expansion_length() + 1, 1);
			res(0) = 0.0;
			for (size_t k = 1; k <= to.get_expansion_length(); ++k)
				res(k) = -std::pow(zy, k-1);
			res *= zydny;
			return -res / two_pi;
		}
	};


	/// \brief P2L operator of the Laplace 2D FMM
	/// \tparam Ny order of derivation
	template <unsigned Ny>
	class p2l
	{
	public:
		typedef cluster_t test_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_2d<double> >,
			0, Ny
		>::trial_input_t trial_input_t;
		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> result_t;

		size_t rows(test_input_t const &to) const
		{
			return to.get_expansion_length() + 1;
		}

		result_t operator()(test_input_t const &to, trial_input_t const &from) const
		{
			return eval(to, from, std::integral_constant<unsigned, Ny>());
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const & from,
			std::integral_constant<unsigned, 0>) const
		{
			using namespace boost::math::double_constants;

			std::complex<double> z = 
				std::complex<double>(from.get_x()(0), from.get_x()(1)) 
				- center2complex(to);

			result_t res(to.get_expansion_length() + 1, 1);
			res(0) = std::log(-z);
			for (size_t k = 1; k <= to.get_expansion_length(); ++k)
				res(k) = -1. / (double(k) * std::pow(z, k));
			return -res / two_pi;
		}

		result_t eval(test_input_t const &to, trial_input_t const & from,
			std::integral_constant<unsigned, 1>) const
		{
			using namespace boost::math::double_constants;
			
			std::complex<double> z =
				std::complex<double>(from.get_x()(0), from.get_x()(1))
				- center2complex(to);

			std::complex<double> zdny(from.get_unit_normal()(0), from.get_unit_normal()(1));

			result_t res(to.get_expansion_length() + 1, 1);
			res(0) = 1. / z;
			for (size_t k = 1; k <= to.get_expansion_length(); ++k)
				res(k) = 1. / std::pow(z, k + 1);
			res *= zdny;
			return -res / two_pi;
		}
	};
	

	/// \brief L2P operator of the Laplace 2D FMM
	/// \tparam Nx order of derivation
	template <unsigned Nx>
	class l2p
	{
	public:
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_2d<double> >,
			Nx, 0
		>::test_input_t test_input_t;

		typedef cluster_t trial_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> complex_result_t;
		typedef complex_result_t result_t;

		size_t cols(trial_input_t const &from) const
		{
			return from.get_expansion_length() + 1;
		}

		result_t operator()(test_input_t const &to, trial_input_t const & from) const
		{
			return eval(to, from, std::integral_constant<unsigned, Nx>());
		}

	private:
		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 0>) const
		{
			std::complex<double> z = std::complex<double>(to.get_x()(0), to.get_x()(1)) - center2complex(from);
			complex_result_t res(1, from.get_expansion_length() + 1);
			for (size_t k = 0; k <= from.get_expansion_length(); ++k)
				res(k) = std::pow(z, k);
			return res;
		}

		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 1>) const
		{
			std::complex<double> z = std::complex<double>(to.get_x()(0), to.get_x()(1)) - center2complex(from);
			std::complex<double> zdnx(to.get_unit_normal()(0), to.get_unit_normal()(1));
			complex_result_t res(1, from.get_expansion_length() + 1);
			for (size_t k = 0; k <= from.get_expansion_length(); ++k)
				res(k) = double(k) * std::pow(z, k-1);
			res *= zdnx;
			return res;
		}
	};
	

	/// \brief M2P operator of the Laplace 2D FMM
	/// \tparam Nx order of derivation
	template <unsigned Nx>
	class m2p
	{
	public:
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_2d<double> >,
			Nx, 0
		>::test_input_t test_input_t;
		typedef cluster_t trial_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> complex_result_t;
		typedef complex_result_t result_t;

		size_t cols(trial_input_t const &from) const
		{
			return from.get_expansion_length() + 1;
		}

		result_t operator()(test_input_t const &to, trial_input_t const &from) const
		{
			return eval(to, from, std::integral_constant<unsigned, Nx>());
		}

	private:
		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 0>) const
		{
			complex_result_t res(1, from.get_expansion_length() + 1);
			std::complex<double> z = std::complex<double>(to.get_x()(0), to.get_x()(1)) - center2complex(from);
			res(0) = std::log(z);
			for (size_t k = 1; k <= from.get_expansion_length(); ++k)
				res(k) = 1. / std::pow(z, k);
			return res;
		}

		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 1>) const
		{
			complex_result_t res(1, from.get_expansion_length() + 1);
			std::complex<double> z = std::complex<double>(to.get_x()(0), to.get_x()(1)) - center2complex(from);
			std::complex<double> zdnx(to.get_unit_normal()(0), to.get_unit_normal()(1));
			res(0) = 1. / z;
			for (size_t k = 1; k <= from.get_expansion_length(); ++k)
				res(k) = -double(k) / std::pow(z, k + 1);
			res *= zdnx;
			return res;
		}
	};


	/// \brief M2M operator of the Laplace 2D FMM
	class m2m
	{
	public:
		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> result_t;
		typedef laplace_2d_fmm::cluster_t cluster_t;

		static size_t unique_idx(cluster_t const &to, cluster_t const &from);

		result_t operator()(cluster_t to, cluster_t from) const;
	};


	/// \brief L2L operator of the Laplace 2D FMM
	class l2l
	{
	public:
		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> result_t;
		typedef laplace_2d_fmm::cluster_t cluster_t;

		static size_t unique_idx(cluster_t const &to, cluster_t const &from);

		result_t operator()(cluster_t to, cluster_t from) const;
	};
	

	/// \brief M2L operator of the Laplace 2D FMM
	class m2l
	{
	public:
		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> result_t;
		typedef laplace_2d_fmm::cluster_t cluster_t;

		static size_t unique_idx(cluster_t const &to, cluster_t const &from);

		result_t operator()(cluster_t to, cluster_t from) const;
	};
	

	/// \brief return an instance of the P2M operator
	/// \tparam Ny order of differentiation w.r.t y
	/// \return an instance of the P2M operator
	template <unsigned Ny>
	p2m<Ny> create_p2m() const
	{
		return p2m<Ny>();
	}
	

	/// \brief return an instance of the P2L operator
	/// \tparam Ny order of differentiation w.r.t y
	/// \return an instance of the P2L operator
	template <unsigned Ny>
	p2l<Ny> create_p2l() const
	{
		return p2l<Ny>();
	}
	

	/// \brief return an instance of the M2P operator
	/// \tparam Nx order of differentiation w.r.t x
	/// \return an instance of the M2P operator
	template <unsigned Nx>
	m2p<Nx> create_m2p() const
	{
		return m2p<Nx>();
	}
	

	/// \brief return an instance of the L2P operator
	/// \tparam Nx order of differentiation w.r.t x
	/// \return an instance of the L2P operator
	template <unsigned Nx>
	l2p<Nx> create_l2p() const
	{
		return l2p<Nx>(); 
	}
	

	/// \brief return an instance of the M2M operator
	/// \return an instance of the M2M operator
	m2m create_m2m() const;
	

	/// \brief return an instance of the L2L operator
	/// \return an instance of the L2L operator
	l2l create_l2l() const;
	

	/// \brief return an instance of the M2L operator
	/// \return an instance of the M2L operator
	m2l create_m2l() const;
};

} // end of namespace fmm
} // namespace NiHu
#endif // LAPLACE_2D_FMM_HPP_INCLUDED
