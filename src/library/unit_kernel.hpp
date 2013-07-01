#ifndef UNIT_KERNEL_HPP_INCLUDED
#define UNIT_KERNEL_HPP_INCLUDED

#include "../bem/kernel.hpp"
#include "../bem/gaussian_quadrature.hpp"

/** \brief the unit kernel returning K(x,y) = 1 for all inputs */
class unit_kernel;

/** \brief traits of the 3D unit kernel */
template<>
struct kernel_traits<unit_kernel>
{
	/** \brief test input type */
	typedef location<space_3d> test_input_t;
	/** \brief trial input type */
	typedef location<space_3d> trial_input_t;
	/** \brief kernel result type */
	typedef double result_t;
	/** \brief quadrature family tag */
	typedef gauss_family_tag quadrature_family_t;
	/** \brief shows if kernel is symmetric */
	static bool const is_symmetric = true;
	/** \brief the singularity order (r^-order) */
	static unsigned const singularity_order = 0;
	/** \brief the singular quadrature order */
	static unsigned const singular_quadrature_order = 0;
};

/**
 * \brief 3D unit kernel
 */
class unit_kernel :
	public kernel_base<unit_kernel>
{
public:
	/** \brief the CRTP base class */
	typedef kernel_base<unit_kernel> base_t;
	
	/**
	 * \brief evaluate kernel at test and trial positions
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \return kernel value K(x,y)
	 */
	result_t operator()(test_input_t const &x, trial_input_t const &y) const
	{
		return 1.0;
	}
	
	/**
	 * \brief evaluate kernel complexity at
	 * \param [in] x test position
	 * \param [in] y trial position
	 * \param [in] s linear element size
	 * \return kernel value K(x,y)
	 */
	unsigned estimate_complexity(test_input_t const &x, trial_input_t const &y, scalar_t const &s) const
	{
		return 0;
	}
};

#endif // UNIT_KERNEL_HPP_INCLUDED

