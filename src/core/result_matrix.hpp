#ifndef RESULT_MATRIX_HPP_INCLUDED
#define RESULT_MATRIX_HPP_INCLUDED

#include "../util/crtp_base.hpp"
#include "../util/couple.hpp"

template <class T>
struct is_result_matrix_impl :
	std::false_type {};

template <class T>
struct is_result_matrix_impl<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > :
	std::true_type {};

template <class Mat>
struct is_result_matrix : is_result_matrix_impl<
	typename std::decay<Mat>::type
> {};

template <class A, class B>
couple<
	typename std::enable_if<is_result_matrix<A>::value, A>::type,
	typename std::enable_if<is_result_matrix<B>::value, B>::type
>
operator,(A &&a, B &&b)
{
	return couple<A, B>(std::forward<A>(a), std::forward<B>(b));
}





#endif // RESULT_MATRIX_HPP_INCLUDED
