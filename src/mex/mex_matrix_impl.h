/**
 * \file mex_matrix_impl.h
 * \brief Implementation of inline functions in mex_matrix.h
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */

#include "mex_matrix.h"

namespace mex {

	//-------------------------------------------------------------------------
	// implementation of  class matrix_base
	//-------------------------------------------------------------------------

	inline matrix_base::matrix_base(size_t rows, size_t cols)
		: m_rows(rows), m_cols(cols)
	{
	}

	inline matrix_base::matrix_base(mxArray const *input)
	{
		m_rows = mxGetM(input);
		m_cols = mxGetN(input);
	}

	inline size_t matrix_base::rows(void) const
	{
		return m_rows;
	}

	inline size_t matrix_base::cols(void) const
	{
		return m_cols;
	}


	//-------------------------------------------------------------------------
	// implementation of  class real_matrix
	//-------------------------------------------------------------------------

	inline real_matrix::real_matrix(size_t rows, size_t cols, mxArray *&output)
		: matrix_base(rows, cols)
	{
		output = mxCreateDoubleMatrix(m_rows, m_cols, mxREAL);
		m_real = mxGetPr(output);
	}

	inline real_matrix::real_matrix(mxArray const *mex_matrix)
		: matrix_base(mex_matrix)
	{
		m_real = mxGetPr(mex_matrix);
	}

	inline double real_matrix::operator() (size_t row, size_t col) const
	{
		return m_real[row + m_rows * col];
	}

	inline double &real_matrix::operator() (size_t row, size_t col)
	{
		return m_real[row + m_rows * col];
	}


	//-------------------------------------------------------------------------
	// implementation of  class complex_matrix
	//-------------------------------------------------------------------------

	inline complex_matrix::complex_matrix(size_t rows, size_t cols, mxArray *&output)
		: matrix_base(rows, cols)
	{
		output = mxCreateDoubleMatrix(m_rows, m_cols, mxCOMPLEX);
		m_real = mxGetPr(output);
		m_imag = mxGetPi(output);
	}

	inline complex_matrix::complex_matrix(mxArray const *mex_matrix)
		: matrix_base(mex_matrix)
	{
		// get real and imaginary data pointers from the Matlab matrix
		m_real = mxGetPr(mex_matrix);
		m_imag = mxGetPi(mex_matrix);
	}

	inline std::complex<double> complex_matrix::operator() (size_t row, size_t col) const
	{
		return std::complex<double>(m_real[row+m_rows*col],
									m_imag[row+m_rows*col]);
	}

	inline index_proxy complex_matrix::operator() (size_t row, size_t col)
	{
		return index_proxy(*this, row, col);
	}


	//-------------------------------------------------------------------------
	// implementation of  class index_proxy
	//-------------------------------------------------------------------------

	inline index_proxy::index_proxy(complex_matrix &matrix, size_t row, size_t col)
		: m_parent(matrix), m_row(row), m_col(col)
	{
	}

	double &index_proxy::real(void) const
	{
		return m_parent.m_real[m_row + m_parent.rows()*m_col];
	}

	double &index_proxy::imag(void) const
	{
		return m_parent.m_imag[m_row + m_parent.rows()*m_col];
	}

	template <class complex_rhs_t>
	inline void index_proxy::operator +=(complex_rhs_t const &data) const
	{
		real() += data.real();
		imag() += data.imag();
	}

} // namespace mex

