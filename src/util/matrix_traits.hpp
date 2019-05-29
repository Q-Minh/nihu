/// \file matrix_traits.hpp
/// \brief compile time properties of matrices
#ifndef MATRIX_TRAITS_HPP_INCLUDED
#define MATRIX_TRAITS_HPP_INCLUDED

#include <complex>

namespace NiHu
{

/** \brief metafunction returning the number of compile time rows
 * \tparam T the matrix type
 * The metafunction returns the number of rows at compile time
 */
template <class T>
struct num_rows
{
	static int const value = T::RowsAtCompileTime;
};

/// \brief specialization of num_rows for the double scalar type
template <>
struct num_rows<double>
{
	static int const value = 1;
};

/// \brief specialization of num_rows for the complex scalar type
template <>
struct num_rows<std::complex<double> >
{
	static int const value = 1;
};


/** \brief metafunction returning the number of compile time columns
 * \tparam T the matrix type
 * The metafunction returns the number of columns at compile time
 */
template <class T>
struct num_cols
{
	static int const value = T::ColsAtCompileTime;
};

/// \brief specialization of num_cols for the double scalar type
template <>
struct num_cols<double>
{
	static int const value = 1;
};

/// \brief specialization of num_cols for the complex scalar type
template <>
struct num_cols<std::complex<double> >
{
	static int const value = 1;
};


/** \brief metafunction returning the scalar type
 * \tparam T the matrix type
 * The metafunction returns the scalar type
 */
template <class T>
struct scalar
{
	typedef typename T::Scalar type;
};


/// \brief specialization of scalar for the double type
template <>
struct scalar<double>
{
	typedef double type;
};


/// \brief specialization of scalar for the complex type
template <>
struct scalar<std::complex<double> >
{
	typedef std::complex<double> type;
};

} // namespace NiHu

#endif // MATRIX_TRAITS_HPP_INCLUDED
