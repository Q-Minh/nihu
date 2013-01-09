#ifndef KERNEL_HPP
#define KERNEL_HPP

#include "bool.hpp"

struct normalKernel;
struct simpleKernel;

template <class C>
struct requires_normal;

template <>
struct requires_normal<normalKernel> : true_ {};

template <>
struct requires_normal<simpleKernel> : false_ {};

#endif

