/**
 * \file includes.h
 * \brief Collection of include files needed by Boonen13
 */

#ifndef INCLUDES_H_INCLUDED
#define INCLUDES_H_INCLUDED

#include <type_traits>

#include <iostream>
#include <complex>
#include <stdexcept>
#include <numeric>

#include <Eigen/Dense>
#include <Eigen/StdVector>
/**
 * \brief macro declaring an Eigen std::vector type with the appropriate allocator
 */
#define EIGENSTDVECTOR(_T) std::vector<_T, Eigen::aligned_allocator<_T> >

#endif // INCLUDES_H_INCLUDED
