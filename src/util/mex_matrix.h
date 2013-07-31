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
 * \details This class provides a common interface to mex matrices.
 * The class stores and can return the matrix size.
 * This class is only used by the dervied classes, and cannot be instantiated.
 */
class matrix_base
{
protected:
	/** \brief output matrix constructor
	 * \details This constructor is used when we create a new matrix in C++ to pass it back to Matlab
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 */
	matrix_base(size_t rows, size_t cols);

	/** \brief input matrix constructor
	 * \details This constructor is used when we receive an input matrix from the Matlab interface
	 * \param [in] input the Matlab pointer to the matrix
	 */
	matrix_base(mxArray const *input);

public:
	/** \brief return number of rows
	 * \return the number of rows
	 */
	size_t rows(void) const;

	/** \brief return number of columns
	 * \return the number of columns
	 */
	size_t cols(void) const;

protected:
	size_t m_rows;		/**< \brief number of rows */
	size_t m_cols;		/**< \brief number of columns */
};



/** \brief container class of a real matrix stored in Matlab format */
class real_matrix : public matrix_base
{
public:
	/** \brief output matrix (allocating) constructor
	 * \details This constructor is used when we create a new matrix in C++ to pass it back to Matlab
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 * \param [out] output pointer to the result matrix  is copied here
	 */
	real_matrix(size_t rows, size_t cols, mxArray *&output);

	/** \brief input matrix constructor
	 * \details This constructor is used when we construct a Matlab-allocated matrix
	 * \param [in] input pointer to the native Matlab matrix format
	 */
	real_matrix(mxArray const *input);

	/** \brief return matrix element
	 * \param [in] row row index
	 * \param [in] col column index
	 * \return matrix element
	 */
	double operator() (size_t row, size_t col) const;

	/** \brief return reference to matrix element
	 * \param [in] row row index
	 * \param [in] col column index
	 * \return reference to matrix element
	 */
	double &operator() (size_t row, size_t col);

protected:
	double *m_real;		/**< \brief array of real data */
};


class index_proxy;

/** \brief container class of a complex matrix stored in Matlab format */
class complex_matrix : public matrix_base
{
public:
	friend class index_proxy;

	/** \brief output matrix (allocating) constructor
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 * \param [out] output pointer to the result matrix  is copied here
	 */
	complex_matrix(size_t rows, size_t cols, mxArray *&output);

	/** \brief input matrix constructor
	 * \param [in] input pointer to the native Matlab matrix format
	 */
	complex_matrix(mxArray const *input);

	/** \brief index operator that returns a complex number
	 * \param row row index
	 * \param col column index
	 * \return complex number
	 */
	std::complex<double> operator() (size_t row, size_t col) const;

	/** \brief index operator that returns a proxy
	 * \param row row index
	 * \param col column index
	 * \return proxy object storing the parent container and the indices
	 */
	index_proxy operator() (size_t row, size_t col);

protected:
	double *m_real;		/**< \brief array of real data */
	double *m_imag;		/**< \brief array of imaginary data */
};


/** \brief index proxy class of a complex matrix */
class index_proxy
{
public:
	/** \brief constructor
	 * \param matrix the parent container class
	 * \param row the row index
	 * \param col the column index
	 */
	index_proxy(complex_matrix &matrix, size_t row, size_t col);

	/** \brief return reference to the real part
	 * \return reference to the real part of the complex element
	 */
	double &real(void) const;

	/** \brief return reference to the imaginary part
	 * \return reference to the imaginary part of the complex element
	 */
	double &imag(void) const;

	/** \brief increment operator
	 * \param data the data to add to the container
	 */
	template <class complex_rhs_t>
	void operator +=(complex_rhs_t const &data) const;

private:
	complex_matrix &m_parent; /**< \brief reference to the parent */
	size_t const m_row; /**< \brief the row index */
	size_t const m_col; /**< \brief the column index */
};

} // namespace mex

#include "mex_matrix_impl.h"

#endif // MEX_MATRIX_H_INCLUDED

