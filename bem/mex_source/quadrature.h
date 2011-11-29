#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "types.h"

void sing_quadr_face_quad(const gauss_t *g,
                          double xi,
                          double eta,
                          int *gnum,
                          double *xiprime,
                          double *etaprime,
                          double *wprime);

void sing_quadr_corner_quad(const gauss_t *g,
                            int corner,
                            int *gnum,
                            double *xiprime,
                            double *etaprime,
                            double *wprime);

void sing_quadr_face_tria(const gauss_t *g,
                          double xi,
                          double eta,
                          int *gnum,
                          double *xiprime,
                          double *etaprime,
                          double *wprime);

void sing_quadr_corner_tria(const gauss_t *g,
                            int corner,
                            int *gnum,
                            double *xiprime,
                            double *etaprime,
                            double *wprime);

#endif
