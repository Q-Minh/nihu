/**
 * \file mex_matrix.hpp
 * \brief A Matlab mex matrix interface
 * \details The interface makes it possible to use Matlab-borne matrices in C++
 * and to pass C++-borne matrices to Matlab. The interface handles real and complex amtrices.
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */
#ifndef MEX_MATRIX_HPP_INCLUDED
#define MEX_MATRIX_HPP_INCLUDED

#include <mex.h>
#include "matrix.h"

#include <complex>

/**
 * \brief base class of a Matlab mex matrix
 * \details This class is only used by the dervied classes, and cannot be instantiated.
 */
class mex_matrix_base
{
protected:
	/** \brief allocating constructor
	 * \details This constructor is used when we create a new matrix in C++ to pass it back to Matlab
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 */
	mex_matrix_base(size_t rows, size_t cols) : m_output(NULL), m_rows(rows), m_cols(cols)
	{
	}

	/** \brief input constructor
	 * \details This constructor is used when we receive an input matrix from the Matlab interface
	 * \param [in] input the Matlab pointer to the matrix
	 */
	mex_matrix_base(mxArray const *input) : m_output(NULL)
	{
		m_rows = mxGetM(input);
		m_cols = mxGetN(input);
	}

public:
	/** \brief return number of rows
	 * \return the number of rows
	 */
	size_t rows(void) const
	{
		return m_rows;
	}

	/** \brief return number of columns
	 * \return the number of columns
	 */
	size_t cols(void) const
	{
		return m_cols;
	}

	/** \brief return the Matlab pointer to the matrix
	 * \details This function is useful when we allocate a matrix in C++ and want to pass its
	 * pointer back to Matlab
	 * \return the Matlab pointer to the matrix
	 */
	mxArray *get_output() const
	{
		return m_output;
	}

protected:
	mxArray *m_output;	/**< \brief pointer to the Matlab matrix */
	size_t m_rows;		/**< \brief number of rows */
	size_t m_cols;		/**< \brief number of columns */
};



/** \brief container class of a real matrix stored in Matlab format */
class mex_real_matrix : public mex_matrix_base
{
public:
	/** \brief output matrix (allocating) constructor
	 * \details This constructor is used when we create a new matrix in C++ to pass it back to Matlab
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 */
	mex_real_matrix(size_t rows, size_t cols) : mex_matrix_base(rows, cols)
	{
		// allocate space for the real matrix
		m_output = mxCreateDoubleMatrix(m_rows, m_cols, mxREAL);
		// get real data pointer from the allocated matrix
		m_real = mxGetPr(m_output);
	}

	/** \brief input matrix constructor
	 * \details This constructor is used when we construct a Matlab-allocated matrix
	 * \param [in] mex_matrix pointer to the native Matlab matrix format
	 */
	mex_real_matrix(mxArray const *mex_matrix) : mex_matrix_base(mex_matrix)
	{
		// get real data pointer from the Matlab matrix
		m_real = mxGetPr(mex_matrix);
	}

	/** \brief return matrix element
	 * \param [in] row row index
	 * \param [in] col column index
	 * \return matrix element
	 */
	double operator() (size_t row, size_t col) const
	{
		return m_real[row + m_rows * col];
	}

	/** \brief return reference to matrix element
	 * \param [in] row row index
	 * \param [in] col column index
	 * \return reference to matrix element
	 */
	double &operator() (size_t row, size_t col)
	{
		return m_real[row + m_rows * col];
	}

protected:
	double *m_real;		/**< \brief array of real data */
};


/** \brief container class of a complex matrix stored in Matlab format */
class mex_complex_matrix : public mex_matrix_base
{
public:
	/** \brief index proxy class of a complex matrix */
	class index_proxy
	{
	public:
		/** \brief constructor
		 * \param matrix the parent container class
		 * \param row the row index
		 * \param col the column index
		 */
		index_proxy(mex_complex_matrix &matrix, size_t row, size_t col)
			: m_matrix(matrix), m_row(row), m_col(col)
		{
		}

		/** \brief increment operator
		 * \param data the data to add to the conatiner
		 */
		void operator +=(std::complex<double> const &data) const
		{
			m_matrix.m_real[m_row+m_matrix.rows()*m_col] += data.real();
			m_matrix.m_imag[m_row+m_matrix.rows()*m_col] += data.imag();
		}

	private:
		mex_complex_matrix &m_matrix; /**< \brief reference to the parent */
		size_t const m_row; /**< \brief the row index */
		size_t const m_col; /**< \brief the column index */
	};

	/** \brief output matrix (allocating) constructor
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 */
	mex_complex_matrix(size_t rows, size_t cols) : mex_matrix_base(rows, cols)
	{
		m_output = mxCreateDoubleMatrix(m_rows, m_cols, mxCOMPLEX);
		m_real = mxGetPr(m_output);
		m_imag = mxGetPi(m_output);
	}

	/** \brief input matrix constructor
	 * \param [in] mex_matrix pointer to the native Matlab matrix format
	 */
	mex_complex_matrix(mxArray const *mex_matrix) : mex_matrix_base(mex_matrix)
	{
		// get real and imaginary data pointers from the Matlab matrix
		m_real = mxGetPr(mex_matrix);
		m_imag = mxGetPi(mex_matrix);
	}

	/** \brief index operator that returns a complex number
	 * \param row row index
	 * \param col column index
	 * \return complex number
	 */
	std::complex<double> ref_index(size_t row, size_t col) const
	{
		return std::complex<double>(m_real[row+m_rows*col],
									m_imag[row+m_rows*col]);
	}

	/** \brief index operator that returns a proxy
	 * \param row row index
	 * \param col column index
	 * \return proxy object storing the parent container and the indices
	 */
	index_proxy ref_index(size_t row, size_t col)
	{
		return index_proxy(*this, row, col);
	}

protected:
	double *m_real;		/**< \brief array of real data */
	double *m_imag;		/**< \brief array of imaginary data */
};

#endif // MEX_MATRIX_HPP_INCLUDED

