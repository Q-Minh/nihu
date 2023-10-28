/** \file laplace_3d_fmm.hpp
 * \brief implementation of the FMM method for the 3D Laplace equation
 */
#ifndef LAPLACE_3D_FMM_HPP_INCLUDED
#define LAPLACE_3D_FMM_HPP_INCLUDED


#include <Eigen/Core>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <complex>

#include "cluster.hpp"
#include "m2l_indices.hpp"
#include "p2p.hpp"

#include "library/laplace_kernel.hpp"
#include "library/laplace_nearly_singular_integrals.hpp"
#include "library/laplace_singular_integrals.hpp"

namespace NiHu
{
namespace fmm
{

/// \brief cluster type of the Laplace 3D FMM
class laplace_3d_cluster;

/// \brief specialisation of cluster traits for the 3D Laplace FMM
template <>
struct cluster_traits<laplace_3d_cluster>
{
private:
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

public:
	/// \brief the space dimension
	static size_t const dimension = 3;
	/// \brief the multipole type
	typedef cvector_t multipole_t;
	/// \brief the local type
	typedef cvector_t local_t;
};

/// \brief cluster type of the Laplace 3D FMM
class laplace_3d_cluster
	: public cluster_base<laplace_3d_cluster>
{
	typedef cluster_base<laplace_3d_cluster> base_t;
	
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

	size_t data_size() const;

	/// \brief return a cleared multipole contribution
	multipole_t zero_multipole() const;

	/// \brief return a cleared local contribution
	local_t zero_local() const;

	Eigen::Index linear_index(Eigen::Index n, Eigen::Index m) const;

private:
	size_t m_expansion_length;
};

/// \brief the fmm for the Laplace equation in 3D
class laplace_3d_fmm
{
public:
	/// \brief the cluster type
	typedef laplace_3d_cluster cluster_t;
	/// \brief the bounding box type
	typedef cluster_t::bounding_box_t bounding_box_t;
	/// \brief the location type
	typedef cluster_t::location_t location_t;

public:
	/// \brief return an instance of the P2P operator
	/// \tparam Nx order of differentiation w.r.t x
	/// \tparam Ny order of differentiation w.r.t y
	/// \return an instance of the P2P operator
	template <int Nx, int Ny>
	fmm::p2p<NiHu::normal_derivative_kernel<
		NiHu::laplace_kernel<NiHu::space_3d<double> >, Nx, Ny
	> >
		create_p2p() const
	{
		typedef NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_3d<double> >, Nx, Ny
		> kernel_t;
		return fmm::p2p<kernel_t>(kernel_t());
	}

	static std::complex<double> Y(size_t n, int m, double theta, double phi)
	{
		using boost::math::spherical_harmonic;
		using namespace boost::math::double_constants;
		
		std::complex<double> res = 
			spherical_harmonic(unsigned(n), std::abs(m), theta, phi) 
			/ std::sqrt((2 * n + 1.0) / (4.0 * pi));
		return (m < 0) ? std::conj(res) : res;
	}

	static double A(int n, int m)
	{
		using boost::math::factorial;
		double res = std::pow(-1, n) / std::sqrt(
			factorial<double>(n - m) * factorial<double>(n + m));
		return res;
	}


	/// \brief P2M operator of the Laplace 3D FMM
	/// \tparam Ny order of derivation
	template <unsigned Ny>
	class p2m
	{
	public:
		typedef cluster_t test_input_t;
		typedef cluster_t::bounding_box_t::location_t x_t;

		typedef typename NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_3d<double> >,
			0, Ny
		>::trial_input_t trial_input_t;

		typedef cluster_t::multipole_t result_t;

		size_t rows(test_input_t const &to) const
		{
			return to.data_size();
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
			
			result_t res(to.data_size(), 1);

			x_t x = from.get_x() - to.get_bounding_box().get_center();
			double r = x.norm();
			double phi = std::atan2(x(1), x(0));
			double theta = std::acos(x(2) / r);

			for (int n = 0; n <= to.get_expansion_length(); ++n)
			{
				double rn = std::pow(r, n);
				for (int m = -n; m <= n; ++m)
				{
					auto idx = to.linear_index(n, m);
					res(idx, 0) = rn * Y(n, -m, theta, phi);
				}
			}

			return res / (4. * pi);
		}

		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 1>) const
		{
			using namespace boost::math::double_constants;
			
			result_t res = to.zero_multipole();
			return res / (4. * pi);
		}
	};


	/// \brief P2L operator of the Laplace 3D FMM
	/// \tparam Ny order of derivation
	template <unsigned Ny>
	class p2l
	{
	public:
		typedef cluster_t::bounding_box_t::location_t x_t;
		typedef cluster_t test_input_t;
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_3d<double> >,
			0, Ny
		>::trial_input_t trial_input_t;
		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> result_t;

		size_t rows(test_input_t const &to) const
		{
			return to.data_size();
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
			
			result_t res(rows(to), 1);

			x_t x = from.get_x() - to.get_bounding_box().get_center();
			double r = x.norm();
			double phi = std::atan2(x(1), x(0));
			double theta = std::acos(x(2) / r);

			for (int n = 0; n <= to.get_expansion_length(); ++n)
			{
				double rn = std::pow(r, -(n + 1));
				for (int m = -n; m <= n; ++m)
				{
					auto idx = to.linear_index(n, m);
					res(idx, 0) = rn * Y(n, -m, theta, phi);
				}
			}

			return res / (4.0 * pi);
		}

		result_t eval(test_input_t const &to, trial_input_t const & from,
			std::integral_constant<unsigned, 1>) const
		{
			throw std::logic_error("Unimplemented");
			return result_t();
		}
	};
	

	/// \brief L2P operator of the Laplace 3D FMM
	/// \tparam Nx order of derivation
	template <unsigned Nx>
	class l2p
	{
	public:
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_3d<double> >,
			Nx, 0
		>::test_input_t test_input_t;
		typedef cluster_t::bounding_box_t::location_t x_t;

		typedef cluster_t trial_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> complex_result_t;
		typedef complex_result_t result_t;

		size_t cols(trial_input_t const &from) const
		{
			return from.data_size();
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
			result_t res(1, from.data_size());

			x_t x = to.get_x() - from.get_bounding_box().get_center();
			double r = x.norm();
			double phi = std::atan2(x(1), x(0));
			double theta = std::acos(x(2) / r);

			for (int n = 0; n <= from.get_expansion_length(); ++n)
			{
				double rn = std::pow(r, n);
				for (int m = -n; m <= n; ++m)
				{
					auto idx = from.linear_index(n, m);
					res(0, idx) = rn * Y(n, m, theta, phi);
				}
			}

			return res;
		}

		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 1>) const
		{
			throw std::logic_error("Unimplemented");
			return result_t();
		}
	};
	

	/// \brief M2P operator of the Laplace 3D FMM
	/// \tparam Nx order of derivation
	template <unsigned Nx>
	class m2p
	{
	public:
		typedef typename NiHu::normal_derivative_kernel<
			NiHu::laplace_kernel<NiHu::space_3d<double> >,
			Nx, 0
		>::test_input_t test_input_t;
		typedef cluster_t trial_input_t;
		typedef Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> complex_result_t;
		typedef complex_result_t result_t;
		typedef cluster_t::bounding_box_t::location_t x_t;

		size_t cols(trial_input_t const &from) const
		{
			return from.data_size();
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
			result_t res(1, from.data_size());

			x_t x = to.get_x() - from.get_bounding_box().get_center();
			double r = x.norm();
			double phi = std::atan2(x(1), x(0));
			double theta = std::acos(x(2) / r);

			for (size_t n = 0; n <= from.get_expansion_length(); ++n)
			{
				double rn = std::pow(r, -int(n+1));
				for (int m = -int(n); m <= int(n); ++m)
				{
					auto idx = from.linear_index(n, m);
					res(0, idx) = rn * Y(n, m, theta, phi);
				}
			}

			return res;
		}

		result_t eval(test_input_t const &to,
			trial_input_t const &from,
			std::integral_constant<unsigned, 1>) const
		{
			throw std::logic_error("Unimplemented");
			return result_t();
		}
	};


	/// \brief M2M operator of the Laplace 3D FMM
	class m2m
	{
	public:
		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> result_t;
		typedef laplace_3d_fmm::cluster_t cluster_t;

		static size_t unique_idx(cluster_t const &to, cluster_t const &from);

		result_t operator()(cluster_t to, cluster_t from) const;
	};


	/// \brief L2L operator of the Laplace 3D FMM
	class l2l
	{
	public:
		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> result_t;
		typedef laplace_3d_fmm::cluster_t cluster_t;

		static size_t unique_idx(cluster_t const &to, cluster_t const &from);

		result_t operator()(cluster_t to, cluster_t from) const;
	};
	

	/// \brief M2L operator of the Laplace 3D FMM
	class m2l
	{
	public:
		typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> result_t;
		typedef laplace_3d_fmm::cluster_t cluster_t;

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
#endif // LAPLACE_3D_FMM_HPP_INCLUDED
