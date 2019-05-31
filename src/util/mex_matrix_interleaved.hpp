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
 * \file mex_matrix_interleaved.hpp
 * \brief A Matlab mex matrix interface
 * \ingroup matlab
 * \details
 * The interface makes it possible to use Matlab-borne matrices in C++ and to
 * create Matlab matrices in C++. The interface handles real and complex
 * matrices in a convenient manner, hiding mex implementation details from the
 * C++ programmer.
 */

#ifndef MEX_MATRIX_INTERLEAVED_HPP_INCLUDED
#define MEX_MATRIX_INTERLEAVED_HPP_INCLUDED

#include "eigen_utils.hpp"
#include "../core/result_matrix.hpp"

#include <complex>
#include <cstddef>
#include <mex.h>
#include <matrix.h>
#include <type_traits>

namespace NiHu
{

/**
 * \brief Matlab mex interface classes
 * \details The mex namespace contains classes and functions that provide an easy-to-use
 * interface to input and output Matlab matrices
 */
namespace mex {

/// \brief metafunction assigning a Matlab class ID to a C type
template <class Scalar>
struct classID;

template <class RealScalar>
struct classID<std::complex<RealScalar> >
	: public classID<RealScalar>
{
};

template <>
struct classID<int>
{
	static mxClassID const value = mxINT32_CLASS;
};

template <>
struct classID<double>
{
	static mxClassID const value = mxDOUBLE_CLASS;
};

template <>
struct classID<float>
{
	static mxClassID const value = mxSINGLE_CLASS;
};



/// \brief metafunction assigning a Matlab complexity to a C type
template <class Scalar>
struct complexity
{
	static mxComplexity const value = mxREAL;
};

template <class RealScalar>
struct complexity<std::complex<RealScalar> >
{
	static mxComplexity const value = mxCOMPLEX;
};

template <class Scalar>
void *get_data_ptr(mxArray const *pa);

template <class Scalar>
Scalar *get_data_pointer(mxArray const *pa)
{
	if (mxIsComplex(pa) != (complexity<Scalar>::value == mxCOMPLEX))
		mexErrMsgIdAndTxt("mex_matrix_interleaved::get_data_pointer",
			"Mex tries to access complex array as real or vica versa");
	return static_cast<Scalar *>(get_data_ptr<Scalar>(pa));
}

template <>
void *get_data_ptr<double>(mxArray const *pa)
{
	return mxGetDoubles(pa);
}

template <>
void *get_data_ptr<float>(mxArray const *pa)
{
	return mxGetSingles(pa);
}

template <>
void *get_data_ptr<int>(mxArray const *pa)
{
	return mxGetInt32s(pa);
}

template <>
void *get_data_ptr<unsigned>(mxArray const *pa)
{
	return mxGetUint32s(pa);
}

template <>
void *get_data_ptr<std::complex<double> >(mxArray const *pa)
{
	return mxGetComplexDoubles(pa);
}

template <>
void *get_data_ptr<std::complex<float> >(mxArray const *pa)
{
	return mxGetComplexSingles(pa);
}

/// \brief Matlab mex matrix
template <class Scalar>
class matrix
	: public Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> >
{
public:
	typedef Scalar scalar_t;
	typedef Eigen::Map<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> > base_t;

	/** \brief output matrix constructor
	 * \details used when a new matrix is created in C++ and passed to Matlab
	 * \param [in] rows number of rows
	 * \param [in] cols number of columns
	 * \param [out] output the mxArray that will be allocated
	 */
	matrix(size_t rows, size_t cols, mxArray*& output)
		: base_t(
			static_cast<scalar_t*>(
				get_data_pointer<scalar_t>(
					output = mxCreateNumericMatrix(
						rows, cols,
						classID<scalar_t>::value,
						complexity<scalar_t>::value
					)
					)
				),
			rows, cols)
	{
	}

	/** \brief input matrix constructor
	 * \details used when an input matrix is received from Matlab
	 * \param [in] input the Matlab pointer to the matrix
	 */
	matrix(mxArray const* input)
		: base_t(
			static_cast<scalar_t*>(get_data_pointer<scalar_t>(input)),
			mxGetM(input),
			mxGetN(input))
	{
	}
};

template <class RealScalar>
using real_matrix = matrix<RealScalar>;

template <class RealScalar>
using complex_matrix = matrix<std::complex<RealScalar> >;

} // end of namespace mex


/** \brief declaring mex::real_matrix as a result matrix type */
template <class RealScalar>
struct is_result_matrix_impl<mex::real_matrix<RealScalar> > : std::true_type {};

/** \brief declaring mex::complex_matrix as a result matrix type */
template <class RealScalar>
struct is_result_matrix_impl<mex::complex_matrix<RealScalar> > : std::true_type {};

} // end of namespace NiHu

#endif /* MEX_MATRIX_INTERLEAVED_HPP_INCLUDED */

