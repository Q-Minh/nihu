/**
* \file elem_traits.hpp
* \brief Basic element type properties
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef ELEM_TRAITS_HPP
#define ELEM_TRAITS_HPP

/**
* \brief Traits class to describe element properties
*/
template <class elemType>
struct elem_traits;

/**
* \brief Three noded Triangle element
*/
class TriaElem;

/**
* \brief Four noded Quad element
*/
class QuadElem;

/**
* \brief Traits class to describe Triangle properties
*/
template<>
struct elem_traits<TriaElem>
{
    enum
    {
        nNodes = 3,
        isLinear = true
    };
};


/**
* \brief Traits class to describe Quad properties
*/
template<>
struct elem_traits<QuadElem>
{
    enum
    {
        nNodes = 4,
        isLinear = false
    };
};

#endif

