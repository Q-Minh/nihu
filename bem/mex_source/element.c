#include "element.h"

void shapefun_quad(double xi, double eta, double *N)
{
    N[0] = (1.0-xi)*(1.0-eta)/4.0;
    N[1] = (1.0+xi)*(1.0-eta)/4.0;
    N[2] = (1.0+xi)*(1.0+eta)/4.0;
    N[3] = (1.0-xi)*(1.0+eta)/4.0;
}

void d_shapefun_quad(double xi, double eta, double *Nxi, double *Neta)
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

void shapefun_tria(double xi, double eta, double *N)
{
    N[0] = (1.0-xi-eta);
    N[1] = xi;
    N[2] = eta;
}

void d_shapefun_tria(double xi, double eta, double *Nxi, double *Neta)
{
    Nxi[0] = -1.0;
    Nxi[1] = 1.0;
    Nxi[2] = 0.0;

    Neta[0] = -1.0;
    Neta[1] = 0.0;
    Neta[2] = 1.0;
}
