#ifndef PRODUCT_TYPE_HPP_INCLUDED
#define PRODUCT_TYPE_HPP_INCLUDED

#include <type_traits>
#include <complex>

#include <Eigen/Dense>

template <class T>
struct is_eigen_matrix
{
	static const bool value = std::is_base_of<Eigen::MatrixBase<T>, T>::value;
};

template <
	class Lhs,
	class Rhs,
	bool isLhsEigen = is_eigen_matrix<Lhs>::value,
	bool isRhsEigen = is_eigen_matrix<Rhs>::value>
struct product_type;

template <class Lhs, class Rhs>
struct product_type<Lhs, Rhs, true, true>
{
	typedef typename Eigen::ProductReturnType<Lhs, Rhs>::Type type;
};

template <class Lhs, class Rhs>
struct product_type<Lhs, Rhs, true, false>
{
	typedef typename Lhs::ScalarMultipleReturnType type;
};

template <class Lhs, class Rhs>
struct product_type<Lhs, std::complex<Rhs>, true, false>
{
	typedef typename Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_multiple2_op<typename Lhs::Scalar, std::complex<Rhs> >,
		const Lhs
	> type;
};

template <class Lhs, class Rhs>
struct product_type<Lhs, Rhs, false, true> : product_type<Rhs, Lhs, true, false> {};


template <>
struct product_type<double, double, false, false>
{
	typedef double type;
};

template <class T>
struct product_type<T, std::complex<T>, false, false>
{
	typedef std::complex<T> type;
};

template <class T>
struct product_type<std::complex<T>, T, false, false>
{
	typedef std::complex<T> type;
};


template <class T, bool isEigen = is_eigen_matrix<T>::value>
struct plain_type
{
	typedef T type;
};

template <class T>
struct plain_type<T, true>
{
	typedef typename T::PlainObject type;
};


#endif // PRODUCT_TYPE_HPP_INCLUDED
