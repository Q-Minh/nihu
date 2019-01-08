/**
 * \file math_constants.hpp 
 * \brief Math constant definitions
 * \ingroup util 
 */

#ifndef MATH_CONSTANTS_HPP_INCLUDED
#define MATH_CONSTANTS_HPP_INCLUDED

#ifdef _MSVC_VER

#define _USING_MATH_DEFINES
#include <cmath>

#else

#define M_PI 3.14159265358979323846

#endif

#define M_EULER_GAMMA 0.57721566490153286060

#endif /* MATH_CONSTANTS_HPP_INCLUDED */
