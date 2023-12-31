/** 
 * @file black_box_fmm.hpp
 * @brief Implementation of the Black Box FMM using Chebyshev interpolation
 * @ingroup bbfmm
 */

#ifndef BLACK_BOX_FMM_HPP_INCLUDED
#define BLACK_BOX_FMM_HPP_INCLUDED

#include "chebyshev_cluster.hpp"
#include "cluster_tree.hpp"
#include "fmm_operator.hpp"
#include "kron_identity.hpp"
#include "m2l_indices.hpp"
#include "nd_cheb.hpp"
#include "p2p.hpp"
#include "fmm_operator.hpp"

#include "nihu/library/location_normal.hpp"
#include "nihu/library/normal_derivative_kernel.hpp"
#include "nihu/util/matrix_traits.hpp"

#include <cstddef>
#include <type_traits>


namespace NiHu
{
namespace fmm
{

template <class Kernel>
struct kernel_derivative_traits
{
	/** @brief Kernel type without normal derivation */
	typedef Kernel kernel_00_t;
	/** @brief Kernel with normal derivative in @i y */
	typedef Kernel kernel_ny_t;
	/** @brief Kernel with normal derivative in @i x */
	typedef Kernel kernel_nx_t;
	static int const nx = 0;
	static int const ny = 0;
};

template <class DistanceDependentKernel, int Nx, int Ny>
struct kernel_derivative_traits<normal_derivative_kernel<DistanceDependentKernel, Nx, Ny> >
{
	typedef normal_derivative_kernel<DistanceDependentKernel, 0, 0> kernel_00_t;
	typedef normal_derivative_kernel<DistanceDependentKernel, 0, Ny> kernel_ny_t;
	typedef normal_derivative_kernel<DistanceDependentKernel, Nx, 0> kernel_nx_t;
	static int const nx = Nx;
	static int const ny = Ny;
};


/** 
 * @brief Black box FMM for a smooth kernel
 * @tparam Kernel Kernel class
 */
template <class Kernel>
class black_box_fmm
{
public:
	/** \brief template parameter as nested type */
	typedef Kernel kernel_t;

	typedef kernel_derivative_traits<kernel_t> derivative_traits_t;

	typedef typename derivative_traits_t::kernel_00_t kernel_00_t;
	typedef typename derivative_traits_t::kernel_nx_t kernel_nx_t;
	typedef typename derivative_traits_t::kernel_ny_t kernel_ny_t;
	static int const Nx = derivative_traits_t::nx;
	static int const Ny = derivative_traits_t::ny;

	/** \brief the space's dimension */
	static size_t const space_dimension = kernel_00_t::space_t::dimension;

	/// \brief the field's dimension
	static size_t const field_dimension = num_rows<typename kernel_00_t::result_t>::value;

	/// \brief the kernel result's scalar type
	typedef typename scalar<typename kernel_00_t::result_t>::type kernel_scalar_t;

	/** \brief the cluster type */
	typedef chebyshev_cluster<space_dimension, kernel_scalar_t, field_dimension> cluster_t;


	/** \brief the m2m operator of the black box fmm */
	class m2m
		: public fmm_operator<m2m_tag>
	{
	public:
		typedef typename black_box_fmm<kernel_t>::cluster_t cluster_t;

		typedef kron_identity<
			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
			field_dimension
		> result_t;

		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return cluster_t::bounding_box_t::dist2idx(
				from.get_bounding_box().get_center(),
				to.get_bounding_box().get_center());
		}

		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			return result_t(chebanterp<double, space_dimension>(
				to.get_chebyshev_order(),
				to.get_bounding_box(),
				from.get_chebyshev_nodes()));
		}
	};


	/** \brief the l2l operator of the black box fmm */
	class l2l
		: public fmm_operator<l2l_tag>
	{
	public:
		typedef typename black_box_fmm<kernel_t>::cluster_t cluster_t;

		typedef kron_identity <
			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
			field_dimension
		> result_t;

		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return cluster_t::bounding_box_t::dist2idx(
				to.get_bounding_box().get_center(),
				from.get_bounding_box().get_center());
		}

		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			return result_t(chebanterp<double, space_dimension>(
				from.get_chebyshev_order(),
				from.get_bounding_box(),
				to.get_chebyshev_nodes()).transpose());
		}
	};


	/** \brief the m2l operator of the black box fmm */
	class m2l
		: public fmm_operator<m2l_tag>
	{
	public:
		typedef typename black_box_fmm<kernel_t>::cluster_t cluster_t;

		typedef Eigen::Matrix<kernel_scalar_t, Eigen::Dynamic, Eigen::Dynamic> result_t;

		m2l(kernel_t const &kernel)
			: m_kernel_00(kernel)
		{
		}

		static size_t unique_idx(cluster_t const &to, cluster_t const &from)
		{
			return m2l_indices<space_dimension>::eval(to.get_bounding_box(), from.get_bounding_box());
		}

		result_t operator()(cluster_t const &to, cluster_t const &from) const
		{
			size_t N = to.get_chebyshev_nodes().cols();
			size_t M = from.get_chebyshev_nodes().cols();
			result_t res(M * field_dimension, N * field_dimension);
			for (size_t i = 0; i < N; ++i)
				for (size_t j = 0; j < M; ++j)
					res.block(i * field_dimension, j * field_dimension, field_dimension, field_dimension) = m_kernel_00(
						typename kernel_t::space_t::location_t(to.get_chebyshev_nodes().col(i)),
						typename kernel_t::space_t::location_t(from.get_chebyshev_nodes().col(j)))
					* Eigen::Matrix<kernel_scalar_t, field_dimension, field_dimension>::Identity();
			return res;
		}

	private:
		kernel_00_t m_kernel_00;
	};


	/** \brief the p2m operator of the black box fmm */
	class p2m
		: public fmm_operator<p2m_tag>
	{
	public:
		typedef cluster_t test_input_t;

		typedef typename std::conditional<
			Ny == 0,
			NiHu::location_input<NiHu::space<double, space_dimension> >,
			NiHu::location_normal_jacobian_input<NiHu::space<double, space_dimension> >
		>::type trial_input_t;

		typedef Eigen::Matrix<double, Eigen::Dynamic, field_dimension> result_t;

		size_t rows(test_input_t const &cluster) const
		{
			return cluster.get_data_size();
		}

		size_t cols(trial_input_t const &) const
		{
			return field_dimension;
		}

		result_t operator()(test_input_t const &to, trial_input_t const &from) const
		{
			return eval(to, from, std::integral_constant<int, Ny>());
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &from, std::integral_constant<int, 0>) const
		{
			Eigen::Matrix<double, Eigen::Dynamic, 1> res0 = chebanterp<double, space_dimension>(
				to.get_chebyshev_order(),
				to.get_bounding_box(), from.get_x());
			result_t res(rows(to), cols(from));
			for (Eigen::Index i = 0; i < res0.rows(); ++i)
				for (Eigen::Index j = 0; j < res0.cols(); ++j)
					res.block(field_dimension * i, field_dimension * j, field_dimension, field_dimension) =
					res0(i, j) * Eigen::Matrix<double, field_dimension, field_dimension>::Identity();
			return res;
		}

		result_t eval(test_input_t const &to, trial_input_t const &from, std::integral_constant<int, 1>) const
		{
			Eigen::Matrix<double, Eigen::Dynamic, 1> res0 = chebanterp_dny<double, space_dimension>(
				to.get_chebyshev_order(),
				to.get_bounding_box(),
				from.get_x(),
				from.get_unit_normal());
			result_t res(rows(to), cols(from));
			for (Eigen::Index i = 0; i < res0.rows(); ++i)
				for (Eigen::Index j = 0; j < res0.cols(); ++j)
					res.block(field_dimension * i, field_dimension * j, field_dimension, field_dimension) =
					res0(i, j) * Eigen::Matrix<double, field_dimension, field_dimension>::Identity();
			return res;
		}
	};


	/** \brief the l2p operator of the black box fmm */
	class l2p
		: public fmm_operator<l2p_tag>
	{
	public:
		typedef typename std::conditional<
			Nx == 0,
			NiHu::location_input<NiHu::space<double, space_dimension> >,
			NiHu::location_normal_jacobian_input<NiHu::space<double, space_dimension> >
		>::type test_input_t;

		typedef cluster_t trial_input_t;

		typedef	Eigen::Matrix<double, field_dimension, Eigen::Dynamic> result_t;

		size_t rows(test_input_t const &) const
		{
			return field_dimension;
		}

		size_t cols(trial_input_t const &cluster) const
		{
			return cluster.get_data_size();
		}

		result_t operator()(test_input_t const &to, trial_input_t const &from) const
		{
			return eval(to, from, std::integral_constant<int, Nx>());
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &from, std::integral_constant<int, 0>) const
		{
			Eigen::Matrix<double, 1, Eigen::Dynamic> res0 = chebanterp<double, space_dimension>(
				from.get_chebyshev_order(),
				from.get_bounding_box(),
				to.get_x()).transpose();

			result_t res(rows(to), cols(from));
			for (Eigen::Index i = 0; i < res0.rows(); ++i)
				for (Eigen::Index j = 0; j < res0.cols(); ++j)
					res.block(field_dimension * i, field_dimension * j, field_dimension, field_dimension) =
					res0(i, j) * Eigen::Matrix<double, field_dimension, field_dimension>::Identity();
			return res;
		}

		result_t eval(test_input_t const &to, trial_input_t const &from, std::integral_constant<int, 1>) const
		{
			Eigen::Matrix<double, 1, Eigen::Dynamic> res0 = chebanterp_dny<double, space_dimension>(
				from.get_chebyshev_order(),
				from.get_bounding_box(),
				to.get_x(),
				to.get_unit_normal()
				).transpose();

			result_t res(rows(to), cols(from));
			for (Eigen::Index i = 0; i < res0.rows(); ++i)
				for (Eigen::Index j = 0; j < res0.cols(); ++j)
					res.block(field_dimension * i, field_dimension * j, field_dimension, field_dimension) =
					res0(i, j) * Eigen::Matrix<double, field_dimension, field_dimension>::Identity();
			return res;
		}
	};


	/** \brief the p2l operator of the black box fmm */
	class p2l
		: public fmm_operator<p2l_tag>
	{
	public:
		typedef cluster_t test_input_t;
		typedef typename kernel_ny_t::trial_input_t trial_input_t;
		typedef Eigen::Matrix<kernel_scalar_t, Eigen::Dynamic, field_dimension> result_t;

		p2l(kernel_t const &kernel)
			: m_kernel_ny(kernel)
		{
		}

		size_t rows(test_input_t const &cluster) const
		{
			return cluster.get_data_size();
		}

		size_t cols(trial_input_t const &) const
		{
			return field_dimension;
		}

		result_t operator()(test_input_t const &to, trial_input_t const &from) const
		{
			return eval(to, from, std::integral_constant<int, Ny>());
		}

	private:
		result_t eval(test_input_t const &to, trial_input_t const &from, std::integral_constant<int, 0>) const
		{
			size_t n = to.get_chebyshev_nodes().cols();
			result_t res(rows(to), cols(from));
			for (size_t i = 0; i < n; ++i)
				res.block(i * field_dimension, 0, field_dimension, field_dimension) =
				m_kernel_ny(to.get_chebyshev_nodes().col(i), from.get_x()) *
				Eigen::Matrix<kernel_scalar_t, field_dimension, field_dimension>::Identity();
			return res;
		}

		result_t eval(test_input_t const &to, trial_input_t const &from, std::integral_constant<int, 1>) const
		{
			size_t n = to.get_chebyshev_nodes().cols();
			result_t res(rows(to), cols(from));
			for (size_t i = 0; i < n; ++i)
				res.block(i * field_dimension, 0, field_dimension, field_dimension) =
				m_kernel_ny(to.get_chebyshev_nodes().col(i), from.get_x(), from.get_unit_normal()) *
				Eigen::Matrix<kernel_scalar_t, field_dimension, field_dimension>::Identity();
			return res;
		}

		kernel_ny_t m_kernel_ny;
	};


	/** \brief the m2p operator of the black box fmm */
	class m2p
		: public fmm_operator<m2p_tag>
	{
	public:
		typedef cluster_t trial_input_t;
		typedef typename kernel_t::test_input_t test_input_t;
		typedef Eigen::Matrix<kernel_scalar_t, field_dimension, Eigen::Dynamic> result_t;

		m2p(kernel_t const &kernel)
			: m_kernel_nx(kernel)
		{
		}

		size_t rows(test_input_t const &) const
		{
			return field_dimension;
		}

		size_t cols(trial_input_t const &cluster) const
		{
			return cluster.get_data_size();
		}

		result_t operator()(test_input_t const &to, trial_input_t const &from) const
		{
			/** \todo the Kronecker identity is computed physilly */
			size_t n = from.get_chebyshev_nodes().cols();
			result_t res(rows(to), cols(from));
			for (size_t i = 0; i < n; ++i)
				res.block(0, i * field_dimension, field_dimension, field_dimension) =
				m_kernel_nx(to.get_x(), from.get_chebyshev_nodes().col(i)) *
				Eigen::Matrix<kernel_scalar_t, field_dimension, field_dimension>::Identity();
			return res;
		}

		/// \todo the derivative is unimplemented here

	private:
		kernel_nx_t m_kernel_nx;
	};


	/** \brief consructor
	 * \param [in] problem the problem instance
	 */
	black_box_fmm(kernel_t const &kernel)
		: m_kernel(kernel)
	{
	}


	/** \brief return the p2p operator instance
	 * \return the p2p operator instance
	 */
	auto create_p2p() const
	{
		return p2p<kernel_t>(m_kernel);
	}

	/** \brief return the p2m operator instance
	 * \return the p2m operator instance
	 */
	auto create_p2m() const
	{
		return p2m();
	}

	/** \brief return the p2l operator instance
	 * \return the p2l operator instance
	 */
	p2l create_p2l() const
	{
		return p2l(m_kernel);
	}

	/** \brief return the m2p operator instance
	 * \return the m2p operator instance
	 */
	auto create_m2p() const
	{
		return m2p(m_kernel);
	}

	/** \brief return the l2p operator instance
	 * \return the l2p operator instance
	 */
	auto create_l2p() const
	{
		return l2p();
	}

	/** \brief return the m2m operator instance
	 * \return the m2m operator instance
	 */
	m2m create_m2m() const
	{
		return m2m();
	}

	/** \brief return the l2l operator instance
	 * \return the l2l operator instance
	 */
	l2l create_l2l() const
	{
		return l2l();
	}

	/** \brief return the m2l operator instance
	 * \return the l2l operator instance
	 */
	m2l create_m2l() const
	{
		return m2l(m_kernel);
	}

private:
	kernel_t m_kernel;
};

} // end of namespace fmm
} // end of namespace NiHu

#endif /* BLACK_BOX_FMM_HPP_INCLUDED */
