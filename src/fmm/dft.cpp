#include "dft.h"
#include <omp.h>

namespace NiHu
{
namespace fmm
{


dft::dft(size_t N, int flag)
	: dft(N, N, flag)
{
}

dft::dft(size_t N, size_t L, int flag)
	: m_N(N)
	, m_L(L)
	, m_flag(flag)
	, m_num_threads(1)
	, m_results(m_num_threads, cvector_t::Zero(m_L,1))
{
	cvector_t in(m_N, 1), out(m_N, 1);
	m_dft_plan = fftw_plan_dft_1d(m_N,
		(fftw_complex *)in.data(),
		(fftw_complex *)out.data(),
		m_flag, FFTW_MEASURE);
}

dft::dft(dft const &other)
	: m_N(other.m_N)
	, m_L(other.m_L)
	, m_flag(other.m_flag)
	, m_num_threads(other.m_num_threads)
	, m_results(other.m_results)
{
	cvector_t in(m_N, 1), out(m_N, 1);
	m_dft_plan = fftw_plan_dft_1d(m_N,
		(fftw_complex *)in.data(),
		(fftw_complex *)out.data(),
		m_flag, FFTW_MEASURE);
}

dft const &dft::operator=(dft const &other)
{
	if (this == &other)
		return *this;

	m_N = other.m_N;
	m_L = other.m_L;
	m_flag = other.m_flag;
	m_num_threads = other.m_num_threads;
	m_results = other.m_results;
	cvector_t in(m_N, 1), out(m_N, 1);
	m_dft_plan = fftw_plan_dft_1d(m_N,
		(fftw_complex *)in.data(),
		(fftw_complex *)out.data(),
		m_flag, FFTW_MEASURE);
	
	return *this;
}

dft::~dft()
{
	fftw_destroy_plan(m_dft_plan);
}

dft::cvector_t const &dft::execute(cvector_t const &input)
{
	size_t idx = omp_get_thread_num();
	fftw_execute_dft(m_dft_plan, (fftw_complex *)input.data(),
		(fftw_complex *)m_results[idx].data());
	return m_results[idx];
}

dft::cvector_t const &dft::get_result() const
{
	size_t idx = omp_get_thread_num();
	return m_results[idx];
}

void dft::set_num_threads(size_t num_threads)
{
	m_num_threads = num_threads;
	m_results.resize(m_num_threads, cvector_t::Zero(m_L, 1));
}

}
} // namespace NiHu

