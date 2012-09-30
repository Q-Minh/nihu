#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include "elem_traits.hpp"
#include "types.h"

double X4[4] = {-1.0, 1.0, 1.0, -1.0};
double Y4[4] = {-1.0, -1.0, 1.0, 1.0};
double X3[3] = {0.0, 1.0, 0.0};
double Y3[3] = {0.0, 0.0, 1.0};

template <class ElemType>
void sing_quadr_face(const gauss_t &g,
                     double xi,
                     double eta,
                     int *gnum,
                     double *xiprime,
                     double *etaprime,
                     double *wprime)
{
    const unsigned nNodes = elem_traits<ElemType>::nNodes;

    double x[4], y[4];
    x[0] = x[1] = xi;
    y[0] = y[1] = eta;

    int k = 0;
    for (int iTria = 0; iTria < nNodes; iTria++)
    {
        if (nNodes == 3)
        {
            x[2] = X3[iTria];
            x[3] = X3[(iTria+1)%nNodes];
            y[2] = Y3[iTria];
            y[3] = Y3[(iTria+1)%nNodes];
        }
        else if (nNodes == 4)
        {
            x[2] = X4[iTria];
            x[3] = X4[(iTria+1)%nNodes];
            y[2] = Y4[iTria];
            y[3] = Y4[(iTria+1)%nNodes];
        }

        for (int iXi = 0; iXi < g.num; iXi++)
        {
            xiprime[k] = etaprime[k] = 0.0;
            double xixi = 0.0, xieta = 0.0, etaxi = 0.0, etaeta = 0.0;
            for (int iNod = 0; iNod < 4; iNod++)
            {
                xiprime[k] += x[iNod] * g.N[iXi+g.num*iNod];
                etaprime[k] += y[iNod] * g.N[iXi+g.num*iNod];

                xixi += x[iNod] * g.Nxi[iXi+g.num*iNod];
                xieta += x[iNod] * g.Neta[iXi+g.num*iNod];
                etaxi += y[iNod] * g.Nxi[iXi+g.num*iNod];
                etaeta += y[iNod] * g.Neta[iXi+g.num*iNod];
            }
            wprime[k] = g.w[iXi] * (xixi*etaeta - xieta*etaxi);

            k++;
        }
    }
    *gnum = k;
}


/*
void sing_quadr_corner_quad(const gauss_t *g, int corner, int *gnum, double *xiprime, double *etaprime, double *wprime)
{
    int k, iXi, iNod, iTria;
    double xixi, xieta, etaxi, etaeta;
    double x[4], y[4];

    x[0] = x[1] = X4[corner];
    y[0] = y[1] = Y4[corner];

    k = 0;
    for (iTria = 0; iTria < 2; iTria++)
    {
        x[2] = X4[(corner+iTria+1)%4];
        x[3] = X4[(corner+iTria+2)%4];
        y[2] = Y4[(corner+iTria+1)%4];
        y[3] = Y4[(corner+iTria+2)%4];

        for (iXi = 0; iXi < g->num; iXi++)
        {
            xiprime[k] = etaprime[k] = 0.0;
            xixi = xieta = etaxi = etaeta = 0.0;
            for (iNod = 0; iNod < 4; iNod++)
            {
                xiprime[k] += x[iNod] * g->N[iXi+g->num*iNod];
                etaprime[k] += y[iNod] * g->N[iXi+g->num*iNod];

                xixi += x[iNod] * g->Nxi[iXi+g->num*iNod];
                xieta += x[iNod] * g->Neta[iXi+g->num*iNod];
                etaxi += y[iNod] * g->Nxi[iXi+g->num*iNod];
                etaeta += y[iNod] * g->Neta[iXi+g->num*iNod];
            }
            wprime[k] = g->w[iXi] * (xixi*etaeta - xieta*etaxi);
            k++;
        }
    }
    *gnum = k;
}

void sing_quadr_corner_tria(const gauss_t *g, int corner, int *gnum, double *xiprime, double *etaprime, double *wprime)
{
    int k, iXi, iNod;
    double xixi, xieta, etaxi, etaeta;
    double x[4], y[4];
    k = 0;

    x[0] = x[1] = X3[corner];
    y[0] = y[1] = Y3[corner];

    x[2] = X3[(corner+1)%3];
    x[3] = X3[(corner+2)%3];
    y[2] = Y3[(corner+1)%3];
    y[3] = Y3[(corner+2)%3];

    for (iXi = 0; iXi < g->num; iXi++)
    {
        xiprime[k] = etaprime[k] = 0.0;
        xixi = xieta = etaxi = etaeta = 0.0;
        for (iNod = 0; iNod < 4; iNod++)
        {
            xiprime[k] += x[iNod] * g->N[iXi+g->num*iNod];
            etaprime[k] += y[iNod] * g->N[iXi+g->num*iNod];

            xixi += x[iNod] * g->Nxi[iXi+g->num*iNod];
            xieta += x[iNod] * g->Neta[iXi+g->num*iNod];
            etaxi += y[iNod] * g->Nxi[iXi+g->num*iNod];
            etaeta += y[iNod] * g->Neta[iXi+g->num*iNod];
        }
        wprime[k] = g->w[iXi] * (xixi*etaeta - xieta*etaxi);
        k++;
    }
    *gnum = k;
}
*/

#endif
