/// \file unit_sphere_interpolator.h
/// \brief implementation of class fmm::interpolator
#ifndef UNIT_SPHERE_INTERPOLATOR_H_INCLUDED
#define UNIT_SPHERE_INTERPOLATOR_H_INCLUDED

#include "unit_sphere.h"
#include <Eigen/Core>
#include <complex>
#include <vector>
#include <fftw3.h>
#include <boost/math/constants/constants.hpp>

namespace NiHu
{
namespace fmm
{

/// \brief class interpolating over the unit sphere
class interpolator
{
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cmatrix_t;
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dmatrix_t;
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> dvector_t;

public:
	interpolator()
		: dft_plan(nullptr)
		, idft_plan(nullptr)
	{
	}

	interpolator(unit_sphere const &Sfrom, unit_sphere const &Sto);

	~interpolator()
	{
		fftw_destroy_plan(dft_plan);
		fftw_destroy_plan(idft_plan);
	}

	interpolator(interpolator const &other)
		: dft_plan(nullptr)
		, idft_plan(nullptr)
	{
		*this = other;
	}

	interpolator(interpolator &&other);

	interpolator const &operator=(interpolator const &other);

	interpolator const &operator=(interpolator &&other);

	cvector_t const &interpolate(cvector_t const &other) const;

private:
	static dmatrix_t Amatrix(int L, int m,
		dvector_t const &theta,
		dvector_t const &xi,
		dvector_t const &wxi);

	static double C(int l, int m);

private:
	size_t m_Lfrom, m_Lto;
	std::vector<dmatrix_t> m_A;
	fftw_plan dft_plan, idft_plan;

	mutable cmatrix_t m_Phi;
	mutable cmatrix_t m_Q;
	mutable cvector_t m_res;
};

} // end of namespace fmm
} // namespace NiHu

#endif
