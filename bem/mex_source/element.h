#ifndef ELEMENT_H
#define ELEMENT_H

void shapefun_quad(double xi,
                   double eta,
                   double *N);

void d_shapefun_quad(double xi,
                     double eta,
                     double *Nxi,
                     double *Neta);
					 
void shapefun_tria(double xi,
                   double eta,
                   double *N);

void d_shapefun_tria(double xi,
                     double eta,
                     double *Nxi,
                     double *Neta);

void shapefun_lin(double xi,
                   double *N);

void d_shapefun_lin(double xi,
                     double *Nxi);

#endif
