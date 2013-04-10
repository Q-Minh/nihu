/**
 * \file product_type.hpp
 * \brief product type calculations
 * \author Peter Fiala fiala@hit.bme.hu and Peter Rucz rucz@hit.bme.hu
 */

#ifndef PRODUCT_TYPE_HPP_INCLUDED
#define PRODUCT_TYPE_HPP_INCLUDED

#include <type_traits>
#include <complex>

#include <Eigen/Dense>

// forward declaration
template <class Derived>
class couple_base;

// forward declaration
template <class First, class Second>
class couple;


/**
 * \brief helper class to decide whether its argument is an eigen matrix or not
 * \details the implementation is not general yet. The class checks if T is a direct
 * derived class of Eigen::MatrixBase. As a second option, the checks if T is Eigen::GeneralProduct
 * \tparam T an arbitrary parameter
 */
template <class T>
struct is_eigen_matrix
{
	static const bool value = std::is_base_of<Eigen::MatrixBase<T>, T>::value;
};

/** \brief specialisation of is_eigen_matrix for the case of an Eigen::GeneralProduct parameter */
template <class Lhs, class Rhs>
struct is_eigen_matrix<Eigen::GeneralProduct<Lhs, Rhs> > : std::true_type {};


/**
 * \brief helper class to decide whether its argument is a couple expression
 * \details The class checks if T is a direct derived class of ::couple_base
 * \tparam T an arbitrary parameter
 */
template <class T>
struct is_couple
{
	static const bool value = std::is_base_of<couple_base<T>, T>::value;
};


/**
 * \brief helper class to determine the product type of two arbitrary types
 * \tparam Lhs left hand side type
 * \tparam Rhs right hand side type
 * \tparam isLhsEigen boolean indicating if Lhs is an eigen matrix
 * \tparam isRhsEigen boolean indicating if Rhs is an eigen matrix
 */
template <class Lhs, class Rhs, bool isLhsEigen, bool isRhsEigen, bool isLhsCouple, bool isRhsCouple>
struct product_type_impl;


template <class Lhs, class Rhs, bool isLhsEigen, bool isRhsEigen>
struct product_type_impl<Lhs, Rhs, isLhsEigen, isRhsEigen, true, false>
{
	typedef couple<
		typename product_type_impl<typename Lhs::first_value_t, Rhs, isLhsEigen, isRhsEigen, false, false>::type,
		typename product_type_impl<typename Lhs::second_value_t, Rhs, isLhsEigen, isRhsEigen, false, false>::type
	> type;
};


template <class Lhs, class Rhs, bool isLhsEigen, bool isRhsEigen>
struct product_type_impl<Lhs, Rhs, isLhsEigen, isRhsEigen, false, true>
{
	typedef couple<
		typename product_type_impl<Lhs, typename Rhs::first_value_t, isLhsEigen, isRhsEigen, false, false>::type,
		typename product_type_impl<Lhs, typename Rhs::second_value_t, isLhsEigen, isRhsEigen, false, false>::type
	> type;
};



/**
 * \brief specialisation of product_type_impl for the case of two eigen matrices
 * \tparam Lhs left hand side type
 * \tparam Rhs right hand side type
 * \details in this case, Eigen::ProductReturnType determines the product type
 */
template <class Lhs, class Rhs>
struct product_type_impl<Lhs, Rhs, true, true, false, false>
{
	typedef typename Eigen::ProductReturnType<Lhs, Rhs>::Type type;
};

/**
 * \brief specialisation of product_type_impl for the case of eigen - non-eigen product
 * \tparam Lhs left hand side type
 * \tparam Rhs right hand side type
 * \details in this case, the result is ScalarMultipleReturnType of the eigen term
 */
template <class Lhs, class Rhs>
struct product_type_impl<Lhs, Rhs, true, false, false, false>
{
	typedef typename Lhs::ScalarMultipleReturnType type;
};


/**
 * \brief specialisation of product_type_impl for the case of eigen - complex non-eigen product
 * \tparam Lhs left hand side type
 * \tparam Rhs right hand side type
 */
template <class Lhs, class Rhs>
struct product_type_impl<Lhs, std::complex<Rhs>, true, false, false, false>
{
	typedef typename Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_multiple2_op<typename Lhs::Scalar, std::complex<Rhs> >,
		const Lhs
	> type;
};


/**
 * \brief specialisation of product_type_impl for the case of non-eigen - eigen product
 * \tparam Lhs left hand side type
 * \tparam Rhs right hand side type
 * \details it is assumed that the product is commutative in this case
 */
template <class Lhs, class Rhs>
struct product_type_impl<Lhs, Rhs, false, true, false, false> : product_type_impl<Rhs, Lhs, true, false, false, false> {};


/** \brief specialisation of product_type_impl for two doubles */
template <>
struct product_type_impl<double, double, false, false, false, false>
{
	typedef double type;
};

/** \brief specialisation of product_type_impl for two double - complex */
template <class T>
struct product_type_impl<T, std::complex<T>, false, false, false, false>
{
	typedef std::complex<T> type;
};

/** \brief specialisation of product_type_impl for two complex - double */
template <class T>
struct product_type_impl<std::complex<T>, T, false, false, false, false>
{
	typedef std::complex<T> type;
};


/**
 * \brief determine product type of two arbitrary types
 * \tparam Lhs left hand side type
 * \tparam Rhs right hand side type
 */
template <class Lhs, class Rhs>
struct product_type : product_type_impl<
	Lhs,
	Rhs,
	is_eigen_matrix<Lhs>::value,
	is_eigen_matrix<Rhs>::value,
	is_couple<Lhs>::value,
	is_couple<Rhs>::value
> {};


#endif // PRODUCT_TYPE_HPP_INCLUDED

