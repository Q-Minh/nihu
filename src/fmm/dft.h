#ifndef FMM_DFT_H_INCLUDED
#define FMM_DFT_H_INCLUDED

#include <complex>
#include <vector>
#include <Eigen/Dense>

#include <fftw3.h>

namespace NiHu
{
namespace fmm
{


class dft
{
public:
	typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> cvector_t;

	dft(size_t N, int flag);

	dft(size_t N, size_t L, int flag);

	dft(dft const &other);

	dft const &operator=(dft const &other);

	~dft();

	cvector_t const &execute(cvector_t const &input);

	cvector_t const &get_result() const;

	void set_num_threads(size_t num_threads);

private:
	size_t m_N, m_L;
	int m_flag;
	size_t m_num_threads;
	std::vector<cvector_t> m_results;
	fftw_plan m_dft_plan;
};

} // end of namespace
} // namespace NiHu

#endif // FMM_DFT_H_INCLUDED
