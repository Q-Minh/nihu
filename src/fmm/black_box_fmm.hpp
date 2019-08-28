/** 
 * @file black_box_fmm.hpp
 * @brief Implementation of the black box FMM using Chebyshev interpolation
 * @ingroup bbfmm
 */

#ifndef BLACK_BOX_FMM_HPP_INCLUDED
#define BLACK_BOX_FMM_HPP_INCLUDED

#include "chebyshev_cluster.hpp"
#include "kron_identity.hpp"
#include "m2l_indices.hpp"
#include "nd_cheb.hpp"
#include "p2p.hpp"

#include "library/location_normal.hpp"
#include "util/matrix_traits.hpp"

#include <type_traits>

namespace NiHu
{
namespace fmm
{

/** 
 * @brief Black box FMM for a smooth kernel
 * @tparam Kernel the kernel class
 */
template <class Problem>
class black_box_fmm
{
public:
	/** \brief template parameter as nested type */
	typedef Problem problem_t;

	/// \brief the base kernel type
	typedef typename problem_t::template kernel<0, 0>::type kernel_00_t;

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
	{
	public:
		typedef typename black_box_fmm<problem_t>::cluster_t cluster_t;

		typedef fmm::kron_identity<
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
			return result_t(fmm::chebanterp<double, space_dimension>(
				to.get_chebyshev_order(),
				to.get_bounding_box(),
				from.get_chebyshev_nodes()));
		}
	};


	/** \brief the l2l operator of the black box fmm */
	class l2l
	{
	public:
		typedef typename black_box_fmm<problem_t>::cluster_t cluster_t;

		typedef fmm::kron_identity <
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
			return result_t(fmm::chebanterp<double, space_dimension>(
				from.get_chebyshev_order(),
				from.get_bounding_box(),
				to.get_chebyshev_nodes()).transpose());
		}
	};


	/** \brief the m2l operator of the black box fmm */
	class m2l
	{
	public:
		typedef typename problem_t::template kernel<0, 0>::type kernel_t;

		typedef typename black_box_fmm<problem_t>::cluster_t cluster_t;

		typedef Eigen::Matrix<kernel_scalar_t, Eigen::Dynamic, Eigen::Dynamic> result_t;

		m2l(kernel_t const &kernel)
			: m_kernel(kernel)
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
					res.block(i * field_dimension, j * field_dimension, field_dimension, field_dimension) = m_kernel(
						typename kernel_t::space_t::location_t(to.get_chebyshev_nodes().col(i)),
						typename kernel_t::space_t::location_t(from.get_chebyshev_nodes().col(j)))
					* Eigen::Matrix<kernel_scalar_t, field_dimension, field_dimension>::Identity();
			return res;
		}

	private:
		kernel_t m_kernel;
	};


	/** \brief the p2m operator of the black box fmm */
	template <unsigned int Ny>
	class p2m
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
			Eigen::Matrix<double, Eigen::Dynamic, 1> res0 = fmm::chebanterp<double, space_dimension>(
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
			Eigen::Matrix<double, Eigen::Dynamic, 1> res0 = fmm::chebanterp_dny<double, space_dimension>(
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
	template <unsigned int Nx>
	class l2p
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
			Eigen::Matrix<double, 1, Eigen::Dynamic> res0 = fmm::chebanterp<double, space_dimension>(
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
			Eigen::Matrix<double, 1, Eigen::Dynamic> res0 = fmm::chebanterp_dny<double, space_dimension>(
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
	template <unsigned int Ny>
	class p2l
	{
	public:
		typedef cluster_t test_input_t;
		typedef typename problem_t::template kernel<0, Ny>::type kernel_t;
		typedef typename kernel_t::trial_input_t trial_input_t;
		typedef Eigen::Matrix<kernel_scalar_t, Eigen::Dynamic, field_dimension> result_t;

		p2l(kernel_t const &kernel)
			: m_kernel(kernel)
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
				m_kernel(to.get_chebyshev_nodes().col(i), from.get_x()) *
				Eigen::Matrix<kernel_scalar_t, field_dimension, field_dimension>::Identity();
			return res;
		}

		result_t eval(test_input_t const &to, trial_input_t const &from, std::integral_constant<int, 1>) const
		{
			size_t n = to.get_chebyshev_nodes().cols();
			result_t res(rows(to), cols(from));
			for (size_t i = 0; i < n; ++i)
				res.block(i * field_dimension, 0, field_dimension, field_dimension) =
				m_kernel(to.get_chebyshev_nodes().col(i), from.get_x(), from.get_unit_normal()) *
				Eigen::Matrix<kernel_scalar_t, field_dimension, field_dimension>::Identity();
			return res;
		}

		kernel_t m_kernel;
	};


	/** \brief the m2p operator of the black box fmm */
	template <unsigned Nx>
	class m2p
	{
	public:
		typedef typename problem_t::template kernel<Nx, 0>::type kernel_t;
		typedef cluster_t trial_input_t;
		typedef typename kernel_t::test_input_t test_input_t;
		typedef Eigen::Matrix<kernel_scalar_t, field_dimension, Eigen::Dynamic> result_t;

		m2p(kernel_t const &kernel)
			: m_kernel(kernel)
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
			size_t n = from.get_chebyshev_nodes().cols();
			result_t res(rows(to), cols(from));
			for (size_t i = 0; i < n; ++i)
				res.block(0, i * field_dimension, field_dimension, field_dimension) =
				m_kernel(to.get_x(), from.get_chebyshev_nodes().col(i)) *
				Eigen::Matrix<kernel_scalar_t, field_dimension, field_dimension>::Identity();
			return res;
		}

	private:
		kernel_t m_kernel;
	};


	/** \brief consructor
	 * \param [in] problem the problem instance
	 */
	black_box_fmm(problem_t const &problem)
		: m_problem(problem)
	{
	}


	/** \brief return the p2p operator instance
	 * \return the p2p operator instance
	 */
	template <unsigned int Nx, unsigned int Ny>
	p2p<typename problem_t::template kernel<Nx, Ny>::type> create_p2p() const
	{
		typedef typename problem_t::template kernel<Nx, Ny>::type kernel_t;
		return p2p<kernel_t>(m_problem.template get_kernel<Nx, Ny>());
	}

	/** \brief return the p2m operator instance
	 * \return the p2m operator instance
	 */
	template <unsigned int Ny = 0>
	p2m<Ny> create_p2m() const
	{
		return p2m<Ny>();
	}

	/** \brief return the p2l operator instance
	 * \return the p2l operator instance
	 */
	template <unsigned int Ny = 0>
	p2l<Ny> create_p2l() const
	{
		return p2l<Ny>(m_problem.template get_kernel<0, Ny>());
	}

	/** \brief return the m2p operator instance
	 * \return the m2p operator instance
	 */
	template <unsigned int Nx = 0>
	m2p<Nx> create_m2p() const
	{
		return m2p<Nx>(m_problem.template get_kernel<Nx, 0>());
	}

	/** \brief return the l2p operator instance
	 * \return the l2p operator instance
	 */
	template <unsigned int Nx = 0>
	l2p<Nx> create_l2p() const
	{
		return l2p<Nx>();
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
		return m2l(m_problem.template get_kernel<0, 0>());
	}

private:
	problem_t m_problem;
};

} // end of namespace fmm
}
#endif // BLACK_BOX_FMM_HPP_INCLUDED
