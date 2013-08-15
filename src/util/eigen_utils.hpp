/**
 * \file eigen_utils.hpp
 * \brief Collection of include files needed by Boonen13
 */

#ifndef EIGEN_UTILS_HPP_INCLUDED
#define EIGEN_UTILS_HPP_INCLUDED

#include <type_traits>
#include <Eigen/Dense>
#include <Eigen/StdVector>

/** \brief metafunction converting a type into an std::vector<T> type with Eigen allocator */
template <class T>
struct EigenStdVector
{
	typedef std::vector<T, Eigen::aligned_allocator<T> > type;
};

/** \brief metafunction determining if its argument is an Eigen expression or not
 * \tparam T the class to investigate
 */
template <class T>
struct is_eigen : std::is_base_of<
	Eigen::EigenBase<typename std::decay<T>::type>,
	typename std::decay<T>::type
> {};

template <class m1, class m2, int t>
struct is_eigen<Eigen::GeneralProduct<m1, m2, t> > : std::true_type {};

#endif // EIGEN_UTILS_HPP_INCLUDED

