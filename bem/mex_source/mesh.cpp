#include "mesh.h"

#include <cmath>
#include "vector.h"

#include "element.hpp"

void init_accelerators(int nnodes,
                       double const nodes[],
                       int nelements,
                       double const elements[],
                       accelerator_t accelerators[])
{
	enum {NDIM = 3, MAXNNODE = 4};

    /* for each element */
    for (int e = 0; e < nelements; e++)
    {
        int elem[MAXNNODE];
        int nvert = (int)elements[e]; /* number of element vertices */
        for (int s = 0; s < nvert; s++)
            elem[s] = (int)elements[e+(s+1)*nelements];

		if (nvert == 3)
		{
        	double nod[3*3];
			for (int j = 0; j < NDIM; j++)
				for (int s = 0; s < 3; s++)
					nod[s+3*j] = nodes[elem[s]+j*nnodes];
			inverse_matrix<TriaElem>(nod, accelerators[e].gradN);
		}

        for (int j = 0; j < NDIM; j++)
        {
            accelerators[e].center[j] = 0.0;
            for (int s = 0; s < nvert; s++)
                accelerators[e].center[j] += nodes[elem[s]+j*nnodes];
            accelerators[e].center[j] /= (double)nvert;
        }

        if (nvert == 3)
        {
            double a[NDIM], b[NDIM], c[NDIM];
            for (int j = 0; j < NDIM; j++)
            {
                a[j] = nodes[elem[1]+j*nnodes]-nodes[elem[0]+j*nnodes];
                b[j] = nodes[elem[2]+j*nnodes]-nodes[elem[0]+j*nnodes];
            }
            cross(a, b, accelerators[e].n0);
        }
        else
        {
            double a[NDIM], b[NDIM], c[NDIM];
            for (int j = 0; j < NDIM; j++)
            {
                a[j] = (nodes[elem[1]+j*nnodes]-nodes[elem[0]+j*nnodes]+nodes[elem[2]+j*nnodes]-nodes[elem[3]+j*nnodes])/4.0;
                b[j] = (nodes[elem[0]+j*nnodes]-nodes[elem[1]+j*nnodes]+nodes[elem[2]+j*nnodes]-nodes[elem[3]+j*nnodes])/4.0;
                c[j] = (nodes[elem[2]+j*nnodes]-nodes[elem[1]+j*nnodes]+nodes[elem[3]+j*nnodes]-nodes[elem[0]+j*nnodes])/4.0;
            }
            cross(a, c, accelerators[e].n0);
            cross(a, b, accelerators[e].nxi);
            cross(b, c, accelerators[e].neta);
        }
    }
}


void init_accelerators2D(int nnodes,
                       double const *nodes,
                       int nelements,
                       double const *elements,
                       accelerator2D_t *accelerators)
{
	enum {NDIM=2, NVERT=2};

    /* for each element */
    for (int e = 0; e < nelements; e++)
    {
        int elem[NVERT];
        for (int s = 0; s < NVERT; s++)
            elem[s] = (int)elements[e+s*nelements];

        for (int j = 0; j < NDIM; j++)
        {
            accelerators[e].center[j] = 0.0;
            for (int s = 0; s < NVERT; s++)
                accelerators[e].center[j] += nodes[elem[s]+j*nnodes];
            accelerators[e].center[j] /= (double)NVERT;
        }

    	double a[NDIM];
		for (int j = 0; j < NDIM; j++)
			a[j] = nodes[elem[1]+j*nnodes]-nodes[elem[0]+j*nnodes];

		accelerators[e].n0[0] = a[1];
		accelerators[e].n0[1] = -a[0];
    }
}

/* ------------------------------------------------------------------------ */
/* Determine Gaussian integration density based on distance between source */
/* and element center */
int gauss_division(double const q[],
                   double const elemcenter[],
                   double const dist[])
{
	enum{NDIM = 3};
    double distance[NDIM];

    for (int j = 0; j < NDIM; j++)
        distance[j] = elemcenter[j] - q[j];
    double d = sqrt(dot(distance,distance));
    if (d < dist[0])
        return 1;
    if (d < dist[1])
        return 2;
    return 3;
}

/* ------------------------------------------------------------------------ */
/* Determine Gaussian integration density based on distance between source */
/* and element center */
int gauss_division2D(double const q[],
                   double const elemcenter[],
                   double const dist[])
{
	enum {NDIM = 2};
    double distance[NDIM];

    for (int j = 0; j < NDIM; j++)
        distance[j] = elemcenter[j] - q[j];
    double d = sqrt(dot(distance,distance));
    if (d < dist[0])
        return 1;
    if (d < dist[1])
        return 2;
    return 3;
}

