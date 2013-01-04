/**
 * \file operator.hpp
 * \brief declaration of operator metafunctions
 * \details This file only declares the binary operators so that they can be overloaded in specific applications
 */
 
#ifndef OPERATOR_HPP
#define OPERATOR_HPP

/** \brief increment operator */
template <class Pos>
struct next;

/** \brief decrement operator */
template <class Pos>
struct prev;

/** \brief binary plus */
template <class Pos1, class Pos2>
struct plus;

/** \brief binary minus */
template <class Pos1, class Pos2>
struct minus;

/** \brief binary multiply */
template <class Pos1, class Pos2>
struct mul;


#endif

