// This file is a part of NiHu, a C++ BEM template library.
//
// Copyright (C) 2012-2014  Peter Fiala <fiala@hit.bme.hu>
// Copyright (C) 2012-2014  Peter Rucz <rucz@hit.bme.hu>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * \file mex_matrix.hpp
 * \brief A Matlab mex matrix interface
 * \ingroup matlab
 * \details The interface makes it possible to use Matlab-borne matrices in C++
 * and to create Matlab matrices in C++. The interface handles real and complex
 * matrices in a convenient manner, hiding mex implementation details from the
 * C++ programmer.
 */
#ifndef MEX_MATRIX_HPP_INCLUDED
#define MEX_MATRIX_HPP_INCLUDED

#include <type_traits>
#include <cstddef>
#include <complex>
#include <mex.h>
#include <matrix.h>

#include "eigen_utils.hpp"
#include "../core/result_matrix.hpp"

/**
 * \brief Matlab mex interface classes
 * \details The mex namespace contains classes and functions that provide an easy-to-use
 * interface to input and output Matlab matrices
 */
namespace mex {

/** \brief base class of a Matlab mex matrix */
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


/** \brief metafunction assigning a Matlab class ID to a C type */
template <class T>
struct classID;

/** \brief specialisation of classID to double */
template <>
struct classID<double>
{
	/** \brief the Matlab class */
	static mxClassID const value = mxDOUBLE_CLASS;
};

/** \brief specialisation of classID to float */
template <>
struct classID<float>
{
	/** \brief the Matlab class */
	static mxClassID const value = mxSINGLE_CLASS;
};


/**
 * brief Container class of a real matrix stored in Matlab format
 * \tparam the real type (float double etc)
 */
template <class T>
class real_matrix :
	public Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
{
public:
	/** \brief the scalar type */
	typedef T scalar_t;
	typedef Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > base_t;

	/** \brief output matrix (allocating) constructor
	 * \details used when a matrix created in C++ is passed to Matlab
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 * \param [out] output pointer to the result matrix  is copied here
	 */
	real_matrix(size_t rows, size_t cols, mxArray *&output)
		: base_t(
			static_cast<scalar_t *>(mxGetData(output = mxCreateNumericMatrix(rows, cols, classID<scalar_t>::value, mxREAL))),
			rows, cols)
	{
	}

	/** \brief input matrix constructor
	 * \details used a Matlab-allocated matrix is read in C++
	 * \param [in] input pointer to the native Matlab matrix format
	 */
	real_matrix(mxArray const *input)
		: base_t(static_cast<scalar_t *>(mxGetData(input)), mxGetM(input), mxGetN(input))
	{
	}
};


/** \brief index proxy class of a complex matrix */
template <class Parent>
class index_proxy
{
public:
	/** \brief the scalar type of the parent class */
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

	/** \brief conversion of element to complex number
	 * \return the complex number
	 */
	operator std::complex<scalar_t>() const
	{
		return std::complex<scalar_t>(real(), imag());
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

	/** \brief assignment operator
	 * \param data the data to assign to the container
	 */
	template <class complex_rhs_t>
	void operator =(complex_rhs_t const &data) const
	{
		real() = data.real();
		imag() = data.imag();
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

	/** \brief increment operator
	 * \param data the data to add to the container
	 */
	void operator +=(scalar_t const &data)
	{
		real() += data;
	}

	/** \brief assignment operator
	 * \param data the data to assign to the container
	 */
	void operator =(scalar_t const &data)
	{
		real() += data;
		imag() = 0.0;
	}

private:
	/** \brief reference to the parent */
	Parent &m_parent;
	/** \brief the row index */
	size_t const m_row;
	/** \brief the column index */
	size_t const m_col;
};


/** \brief Container class of a complex matrix stored in Matlab format */
template <class T>
class complex_matrix :
	public matrix_base
{
public:
	friend class index_proxy<complex_matrix<T> >;
	/** \brief the real scalar type */
	typedef T scalar_t;

	/** \brief output matrix (allocating) constructor
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 * \param [out] output pointer to the result matrix  is copied here
	 */
	complex_matrix(size_t rows, size_t cols, mxArray *&output)
		: matrix_base(rows, cols)
	{
		output = mxCreateNumericMatrix(
			m_rows, m_cols, classID<scalar_t>::value, mxCOMPLEX);
		m_real = static_cast<scalar_t *>(mxGetData(output));
		m_imag = static_cast<scalar_t *>(mxGetImagData(output));
	}

	/** \brief input matrix constructor
	 * \param [in] input pointer to the native Matlab matrix format
	 */
	complex_matrix(mxArray const *input)
		: matrix_base(input)
	{
		// get real and imaginary data pointers from the Matlab matrix
		m_real = static_cast<scalar_t *>(mxGetData(input));
		m_imag = static_cast<scalar_t *>(mxGetImagData(input));
	}

	/** \brief index operator that returns a complex number
	 * \param row row index
	 * \param col column index
	 * \return complex number
	 */
	std::complex<scalar_t> operator() (size_t row, size_t col) const
	{
		return std::complex<scalar_t>(
			m_real[row+m_rows*col], m_imag[row+m_rows*col]);
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


/** \brief declaring mex::real_matrix as a result matrix type */
template <class T>
struct is_result_matrix_impl<mex::real_matrix<T> > : std::true_type {};

/** \brief declaring mex::complex_matrix as a result matrix type */
template <class T>
struct is_result_matrix_impl<mex::complex_matrix<T> > : std::true_type {};

#endif // MEX_MATRIX_HPP_INCLUDED

