/** \file helmholtz_3d_hf_level_data.h
 * \brief level data of the the helmholtz 3d high frequency fmm
 */
#ifndef FMM_HELMHOLTZ_3D_HF_LEVEL_DATA_H_INCLUDED
#define FMM_HELMHOLTZ_3D_HF_LEVEL_DATA_H_INCLUDED

#include <Eigen/Dense>

#include "unit_sphere.h"
#include "unit_sphere_interpolator.h"

#include <vector>

namespace NiHu
{
namespace fmm
{

/** \brief level data of the helmholtz 3d hf fmm */
class helmholtz_3d_hf_level_data
{
public:
	/** \brief complex dynamic Eigen vector */
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

	/** \brief constructor
	 * Sets vectors to size 1, as they need to be expanded based on first element
	 */
	helmholtz_3d_hf_level_data()
		: m_interp_ups(1)
		, m_interp_dns(1)
	{
	}

	/** \brief set the expansion length
	 * \param [in] L the expansion  length
	 */
	void set_expansion_length(size_t L);

	/** \brief return the expansion length
	 * \return the expansion length
	 */
	size_t get_expansion_length() const;

	/** \brief return the stored unit sphere
	 * \return the stored unit sphere
	 */
	unit_sphere const &get_unit_sphere() const;

	/** \brief set the number of threads
	 * \param [in] n  number of threads
	 */
	void set_num_threads(size_t n);

	/** \brief interpolate a function on the unit sphere up to this level
	 * \param [in] x the function to be interpolated
	 * \return the interpolated function
	 */
	cvector_t const &interp_up(cvector_t const &x) const;

	/** \brief interpolate a function on the unit sphere down to this level
	 * \param [in] x the function to be interpolated
	 * \return the interpolated function
	 */
	cvector_t const &interp_down(cvector_t const &x) const;

	/** \brief set the up-interpolator object (for all threads)
	 * \param [in] interp the interpolator
	 */
	void set_interp_up(interpolator const &interp);

	/** \brief set the down-interpolator object (for all threads)
	 * \param [in] interp the interpolator
	 */
	void set_interp_dn(interpolator const &interp);

private:
	size_t m_expansion_length; // L
	unit_sphere m_unit_sphere;
	std::vector<interpolator> m_interp_ups;
	std::vector<interpolator> m_interp_dns;
};

} // namespace fmm
} // namespace NiHu

#endif // FMM_HELMHOLTZ_3D_HF_LEVEL_DATA_H_INCLUDED
