/**
* \file element.hpp
* \brief Element shape functions
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "elem_traits.hpp"
#include "vector.h"

/**
* \brief Template function to return shape functions of an element
*/
template <class ElemType>
void shapefun(double xi, double eta, double N[]);

/**
* \brief Template function to return derivatives of shape functions of an element
*/
template <class ElemType>
void d_shapefun(double xi, double eta, double Nxi[], double Neta[]);

/**
* \brief Shape functions of a Quad element
* \param xi first coordinate on base element
* \param eta second coordinate on base element
* \param N shape function vector
*/
template <>
inline void shapefun<QuadElem>(double xi, double eta, double N[])
{
    N[0] = (1.0-xi)*(1.0-eta)/4.0;
    N[1] = (1.0+xi)*(1.0-eta)/4.0;
    N[2] = (1.0+xi)*(1.0+eta)/4.0;
    N[3] = (1.0-xi)*(1.0+eta)/4.0;
}

/**
* \brief Shape functions of a Tria element
* \param xi first coordinate on base element
* \param eta second coordinate on base element
* \param N shape function vector
*/
template <>
inline void shapefun<TriaElem>(double xi, double eta, double N[])
{
    N[0] = (1.0-xi-eta);
    N[1] = xi;
    N[2] = eta;
}

/**
* \brief Shape function derivatives of a Quad element
* \param xi first coordinate on base element
* \param eta second coordinate on base element
* \param Nxi shape function derivatives with respect to xi
* \param Neta shape function derivatives with respect to eta
*/
template <>
inline void d_shapefun<QuadElem>(double xi, double eta, double Nxi[], double Neta[])
{
    Nxi[0] = (-1.0)*(1.0-eta)/4.0;
    Nxi[1] = (+1.0)*(1.0-eta)/4.0;
    Nxi[2] = (+1.0)*(1.0+eta)/4.0;
    Nxi[3] = (-1.0)*(1.0+eta)/4.0;

    Neta[0] = (1.0-xi)*(-1.0)/4.0;
    Neta[1] = (1.0+xi)*(-1.0)/4.0;
    Neta[2] = (1.0+xi)*(+1.0)/4.0;
    Neta[3] = (1.0-xi)*(+1.0)/4.0;
}

/**
* \brief Shape function derivatives of a Tria element
* \param xi first coordinate on base element
* \param eta second coordinate on base element
* \param Nxi shape function derivatives with respect to xi
* \param Neta shape function derivatives with respect to eta
*/
template <>
inline void d_shapefun<TriaElem>(double xi, double eta, double Nxi[], double Neta[])
{
    Nxi[0] = -1.0;
    Nxi[1] = 1.0;
    Nxi[2] = 0.0;

    Neta[0] = -1.0;
    Neta[1] = 0.0;
    Neta[2] = 1.0;
}


template <class ElemType>
void inverse_matrix(double *nodes, double *gradN, double xi = 0.0, double eta = 0.0)
{
    enum {NDIM = 3};
    const bool isLinear = elem_traits<ElemType>::isLinear;
    const unsigned nNodes = elem_traits<ElemType>::nNodes;

	double dNxi[nNodes], dNeta[nNodes];
	double drxi[NDIM], dreta[NDIM];

	d_shapefun<ElemType>(xi, eta, dNxi, dNeta);

	for (int j = 0; j < NDIM; j++)
	{
		drxi[j] = dreta[j] = 0.0;
		for (int s = 0; s < nNodes; s++)
		{
			drxi[j] += nodes[s+nNodes*j] * dNxi[s];
			dreta[j] += nodes[s+nNodes*j] * dNeta[s];
		}
	}

	double norm[NDIM], c[NDIM], d[NDIM];
	cross(drxi, dreta, norm);
	cross(dreta, norm, c);
	cross(norm, drxi, d);

	double b = dot(norm, norm);

	for (int j = 0; j < NDIM; ++j)
		for (int s = 0; s < nNodes; s++)
			gradN[j+NDIM*s] = (c[j] * dNxi[s] + d[j] * dNeta[s]) / b;
}

#endif

