/** \file helmholtz_2d_wb_cluster.h
 * \brief declaration of class helmholtz_2d_wb_cluster
 */

#ifndef HELMHOLTZ_2D_WB_CLUSTER_H_INCLUDED
#define HELMHOLTZ_2D_WB_CLUSTER_H_INCLUDED

#include "cluster.hpp"
#include "helmholtz_2d_wb_level_data.h"

#include <Eigen/Dense>

namespace NiHu
{
namespace fmm
{

/** \brief the cluster type of the helmholtz 2d wide band fmm */
class helmholtz_2d_wb_cluster;

/// \brief specialisation of fmm::cluster_traits for the 2d Helmholtz wb fmm
template <>
struct cluster_traits<helmholtz_2d_wb_cluster>
{
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;
	static size_t const dimension = 2;
	typedef cvector_t multipole_t;
	typedef cvector_t local_t;
};


/** \brief the cluster type of the helmholtz 2d wide band fmm */
class helmholtz_2d_wb_cluster
	: public cluster_base<helmholtz_2d_wb_cluster>
{
public:
	typedef cluster_base<helmholtz_2d_wb_cluster> base_t;
	typedef typename base_t::multipole_t multipole_t;
	typedef typename base_t::local_t local_t;

	/** \brief set the pointer to the level data
	 * \param [in] p pointer to the level data
	 */
	void set_p_level_data(helmholtz_2d_wb_level_data const *p);

	/** \brief return the level data
	 * \return the level data
	 */
	helmholtz_2d_wb_level_data const &get_level_data() const;

	/** \brief returne a zero multipole
	 * \return zero multipole
	 */
	multipole_t zero_multipole() const;

	/** \brief returne a zero local
	 * \return zero local
	 */
	local_t zero_local() const;

private:
	helmholtz_2d_wb_level_data const *m_p_level_data;
};

} // end of namespace fmm
} // namespace NiHu

#endif //  HELMHOLTZ_2D_WB_CLUSTER_H_INCLUDED
