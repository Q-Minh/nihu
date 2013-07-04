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
#include <utility> // pair

#include <Eigen/Dense>
#include <Eigen/StdVector>

template <class T>
struct EigenStdVector
{
	typedef std::vector<T, Eigen::aligned_allocator<T> > type;
};


namespace formalism
{
	struct collocational {};
	struct general {};
	struct full_dirac {};
}


#endif // INCLUDES_H_INCLUDED
