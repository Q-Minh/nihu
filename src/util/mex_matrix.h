/**
 * \file mex_matrix.h
 * \brief A Matlab mex matrix interface
 * \details The interface makes it possible to use Matlab-borne matrices in C++
 * and to create Matlab matrices in C++. The interface handles real and complex matrices
 * in a convenient manner, hiding mex implementation details from the C++ programmer.
 * \author Peter Fiala fiala@hit.bme.hu, Peter Rucz rucz@hit.bme.hu
 */
#ifndef MEX_MATRIX_H_INCLUDED
#define MEX_MATRIX_H_INCLUDED

#include <cstddef>
#include <complex>
#include <mex.h>
#include <matrix.h>

/**
 * \brief Matlab mex interface classes
 * \details The mex namespace contains classes and functions that provide an easy-to-use
 * interface to input and output Matlab matrices
 */
namespace mex {

/**
 * \brief base class of a Matlab mex matrix
 */
class matrix_base
{
protected:
	/** \brief output matrix constructor
	 * \details used when a new matrix is created in C++ and passed to Matlab
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 */
	matrix_base(size_t rows, size_t cols)
		: m_rows(rows), m_cols(cols)
	{
	}

	/** \brief input matrix constructor
	 * \details used when an input matrix is received from Matlab
	 * \param [in] input the Matlab pointer to the matrix
	 */
	matrix_base(mxArray const *input)
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

protected:
	/** \brief number of rows */
	size_t m_rows;
	/** \brief number of columns */
	size_t m_cols;
};



/** \brief container class of a real matrix stored in Matlab format */
template <class T>
class real_matrix :
	public matrix_base
{
public:
	typedef T scalar_t;

	/** \brief output matrix (allocating) constructor
	 * \details used when a matrix created in C++ is passed to Matlab
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 * \param [out] output pointer to the result matrix  is copied here
	 */
	real_matrix(size_t rows, size_t cols, mxArray *&output)
		: matrix_base(rows, cols)
	{
		output = mxCreateDoubleMatrix(m_rows, m_cols, mxREAL);
		m_real = mxGetPr(output);
	}

	/** \brief input matrix constructor
	 * \details used a Matlab-allocated matrix is read in C++
	 * \param [in] input pointer to the native Matlab matrix format
	 */
	real_matrix(mxArray const *input)
		: matrix_base(input)
	{
		m_real = mxGetPr(input);
	}

	/** \brief return matrix element
	 * \param [in] row row index
	 * \param [in] col column index
	 * \return matrix element
	 */
	scalar_t operator() (size_t row, size_t col) const
	{
		return m_real[row + m_rows * col];
	}

	/** \brief return reference to matrix element
	 * \param [in] row row index
	 * \param [in] col column index
	 * \return reference to matrix element
	 */
	scalar_t &operator() (size_t row, size_t col)
	{
		return m_real[row + m_rows * col];
	}

protected:
	/** \brief array of real data */
	scalar_t *m_real;
};


/** \brief index proxy class of a complex matrix */
template <class Parent>
class index_proxy
{
public:
	typedef typename Parent::scalar_t scalar_t;

	/** \brief constructor
	 * \param matrix the parent container class
	 * \param row the row index
	 * \param col the column index
	 */
	index_proxy(Parent &matrix, size_t row, size_t col)
		: m_parent(matrix), m_row(row), m_col(col)
	{
	}

	/** \brief return reference to the real part
	 * \return reference to the real part of the complex element
	 */
	scalar_t &real(void) const
	{
		return m_parent.m_real[m_row + m_parent.rows()*m_col];
	}

	/** \brief return reference to the imaginary part
	 * \return reference to the imaginary part of the complex element
	 */
	scalar_t &imag(void) const
	{
		return m_parent.m_imag[m_row + m_parent.rows()*m_col];
	}

	/** \brief increment operator
	 * \param data the data to add to the container
	 */
	template <class complex_rhs_t>
	void operator +=(complex_rhs_t const &data) const
	{
		real() += data.real();
		imag() += data.imag();
	}

private:
	Parent &m_parent; /**< \brief reference to the parent */
	size_t const m_row; /**< \brief the row index */
	size_t const m_col; /**< \brief the column index */
};


/** \brief container class of a complex matrix stored in Matlab format */
template <class T>
class complex_matrix :
	public matrix_base
{
public:
	typedef T scalar_t;

	/** \brief output matrix (allocating) constructor
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 * \param [out] output pointer to the result matrix  is copied here
	 */
	complex_matrix(size_t rows, size_t cols, mxArray *&output)
		: matrix_base(rows, cols)
	{
		output = mxCreateDoubleMatrix(m_rows, m_cols, mxCOMPLEX);
		m_real = mxGetPr(output);
		m_imag = mxGetPi(output);
	}

	/** \brief input matrix constructor
	 * \param [in] input pointer to the native Matlab matrix format
	 */
	complex_matrix(mxArray const *input)
		: matrix_base(input)
	{
		// get real and imaginary data pointers from the Matlab matrix
		m_real = mxGetPr(input);
		m_imag = mxGetPi(input);
	}

	/** \brief index operator that returns a complex number
	 * \param row row index
	 * \param col column index
	 * \return complex number
	 */
	std::complex<scalar_t> operator() (size_t row, size_t col) const
	{
		return std::complex<scalar_t>(m_real[row+m_rows*col],
									m_imag[row+m_rows*col]);
	}

	/** \brief index operator that returns a proxy
	 * \param row row index
	 * \param col column index
	 * \return proxy object storing the parent container and the indices
	 */
	index_proxy<complex_matrix> operator() (size_t row, size_t col)
	{
		return index_proxy<complex_matrix>(*this, row, col);
	}

protected:
	/** \brief array of real data */
	scalar_t *m_real;
	/** \brief array of imaginary data */
	scalar_t *m_imag;
};


} // namespace mex


#include "../bem/result_matrix.hpp"

template <class T>
struct is_result_matrix_impl<mex::real_matrix<T> > : std::true_type {};

template <class T>
struct is_result_matrix_impl<mex::complex_matrix<T> > : std::true_type {};

#endif // MEX_MATRIX_H_INCLUDED

