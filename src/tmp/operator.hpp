/**
 * \file operator.hpp
 * \brief declaration of operator metafunctions
 * \details This file only declares the binary operators so that they can be overloaded in specific applications
 */
 
#ifndef OPERATOR_HPP
#define OPERATOR_HPP

/** \brief increment operator */
template <class A>
struct next;

/** \brief decrement operator */
template <class A>
struct prev;

/** \brief binary plus */
template <class A, class B>
struct plus;

/** \brief binary minus */
template <class A, class B>
struct minus;

/** \brief binary multiply */
template <class A, class B>
struct mul;

/** \brief binary less than */
template <class A, class B>
struct less;

#endif

