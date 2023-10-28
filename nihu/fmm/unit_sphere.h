/** \file unit_sphere.h
 * \brief declaration of class fmm::unit_sphere
 */
#ifndef UNIT_SPHERE_HPP_INCLUDED
#define UNIT_SPHERE_HPP_INCLUDED

#include "Eigen/Dense"

namespace NiHu
{
namespace fmm
{

/** \brief class performing interpolation and integration over the unit sphere */
class unit_sphere
{
	/// \brief complex column vector type
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;
	/// \brief real column vector type
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dvector_t;
	/// \brief complex matrix type
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cmatrix_t;
	/// \brief real matrix type
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dmatrix_t;

	size_t m_P;	// quadrature density
	dvector_t theta, phi, w, wtheta;
	Eigen::Matrix<double, 3, Eigen::Dynamic> s;

public:
	unit_sphere()
	{
	}

	/** \brief constructor from quadrature density
	 * \param [in] P quadrature density (number of quadrature points in theta direction)
	 */
	unit_sphere(size_t P);

	/** \brief return quadrature density
	 * \return quadrature density
	 */
	size_t get_P(void) const
	{
		return m_P;
	}

	/** \brief return quadrature weights
	 * \return quadrature weights
	 */
	dvector_t const &get_w(void) const
	{
		return w;
	}

	/** \brief return theta quadrature points
	 * \return theta quadrature points
	 */
	dvector_t const &get_theta(void) const
	{
		return this->theta;
	}

	/** \brief return theta quadrature weights
	 * \return theta quadrature weights
	 */
	dvector_t const &get_wtheta(void) const
	{
		return this->wtheta;
	}

	/** \brief return phi quadrature points
	 * \return phi quadrature points
	 */
	dvector_t const &get_phi(void) const
	{
		return phi;
	}

	/** \brief return Cartesian quadrature points
	 * \return Cartesian quadrature points
	 */
	Eigen::Matrix<double, 3, Eigen::Dynamic> const &get_s(void) const
	{
		return s;
	}

	/** \brief perform spherical transform
	 * \param [in] f vector of function samples
	 * \param [in] L expansion length
	 * \return spherical harmonc expansion coefficients
	 */
	cmatrix_t spht(cvector_t const &f, int L) const;

	/** \brief perform inverse spherical transform
	 * \param [in] Flm spherical harmonc expansion coefficients
	 * \return vector of function samples
	 */
	cvector_t ispht(cmatrix_t const &Flm) const;
};

} // namespace fmm
} // namespace NiHu

#endif // UNIT_SPHERE_HPP_INCLUDED
