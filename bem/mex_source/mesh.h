/**
* \file mesh.h
* \brief Header file of mesh.cpp
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef MESH_H
#define MESH_H

#include "types.h"

void init_accelerators(int nnodes,
                       double const nodes[],
                       int nelements,
                       double const elements[],
                       accelerator_t accelerators[]);

void init_accelerators2D(int nnodes,
                         double const nodes[],
                         int nelements,
                         double const elements[],
                         accelerator2D_t accelerators[]);

int gauss_division(double const q[],
                   double const elemcenter[],
                   double const dist[]);

int gauss_division2D(double const q[],
                     double const elemcenter[],
                     double const dist[]);

#endif

