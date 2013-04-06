#ifndef PRODUCT_TYPE_HPP_INCLUDED
#define PRODUCT_TYPE_HPP_INCLUDED

#include <complex>
#include <Eigen/Dense>
typedef std::complex<double> dcomplex;

template <class Lhs, class Rhs>
struct product_type;

template <>
struct product_type<double, double>
{
	typedef double type;
};

template <class T>
struct product_type<T, std::complex<T> >
{
	typedef std::complex<T> type;
};

template <class T>
struct product_type<std::complex<T>, T>
{
	typedef std::complex<T> type;
};

template <class T, int Rows, int Cols>
struct product_type<Eigen::Matrix<T, Rows, Cols>, std::complex<T> >
{
	typedef Eigen::Matrix<std::complex<double>, Rows, Cols> type;
};

template <class T, int Rows, int Cols>
struct product_type<std::complex<T>, Eigen::Matrix<T, Rows, Cols> >
{
	typedef Eigen::Matrix<std::complex<double>, Rows, Cols> type;
};

template <class T1, int Rows, int Common, class T2, int Cols>
struct product_type<Eigen::Matrix<T1, Rows, Common>, Eigen::Matrix<T2, Common, Cols> >
{
	typedef Eigen::Matrix<typename product_type<T1, T2>::type, Rows, Cols> type;
};



#endif // PRODUCT_TYPE_HPP_INCLUDED
