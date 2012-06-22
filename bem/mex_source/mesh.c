#include "mesh.h"
#include <math.h>
#include "vector.h"

void inverse_matrix_tri(double *nodes, double *gradN)
{
#define NVERT 3
#define NDIM 3
	double dNxi[NVERT] = {-1.0, 1.0, 0.0};
	double dNeta[NVERT] = {-1.0, 0.0, 1.0};
	double drxi[NDIM], dreta[NDIM], norm[NDIM], b;
	double c[NDIM], d[NDIM];
	int j, s;
	
	for (j = 0; j < NDIM; j++)
	{
		drxi[j] = dreta[j] = 0.0;
		for (s = 0; s < NVERT; s++)
		{
			drxi[j] += nodes[s+NVERT*j] * dNxi[s];
			dreta[j] += nodes[s+NVERT*j] * dNeta[s];
		}
	}
	
	cross(drxi, dreta, norm);
	cross(dreta, norm, c);
	cross(norm, drxi, d);
	
	b = dot(norm, norm);
	
	for (j = 0; j < NDIM; ++j) 
		for (s = 0; s < NVERT; s++)
			gradN[j+NDIM*s] = (c[j] * dNxi[s] + d[j] * dNeta[s]) / b;
#undef NVERT
#undef NDIM
}

void inverse_matrix_quad(const double *nodes, double xi, double eta, double *gradN)
{
#define NVERT 4
#define NDIM 3
	double dNxi[NVERT], dNeta[NVERT];				/* Shape functions' derivatives */
	double drxi[NDIM], dreta[NDIM], norm[NDIM], b;
	double c[NDIM], d[NDIM];
	int j, s;
	
	/* Evaluate shape functions */
	d_shapefun_quad(xi, eta, dNxi, dNeta);
	
	for (j = 0; j < NDIM; j++)
	{
		drxi[j] = dreta[j] = 0.0;
		for (s = 0; s < NVERT; s++)
		{
			drxi[j] += nodes[s+NVERT*j] * dNxi[s];
			dreta[j] += nodes[s+NVERT*j] * dNeta[s];
		}
	}
	
	cross(drxi, dreta, norm);
	cross(dreta, norm, c);
	cross(norm, drxi, d);

	b = dot(norm, norm);
	
	for (j = 0; j < NDIM; ++j)
		for (s = 0; s < NVERT; s++)
			gradN[j+NDIM*s] = (c[j] * dNxi[s] + d[j] * dNeta[s]) / b;
#undef NVERT
#undef NDIM
}

void init_accelerators(int nnodes,
                       const double *nodes,
                       int nelements,
                       const double *elements,
                       accelerator_t *accelerators)
{
    int e, j, s, nvert;
    int elem[4];
    double a[3], b[3], c[3];
	double nod[3*3];

    /* for each element */
    for (e = 0; e < nelements; e++)
    {
        nvert = (int)elements[e]; /* number of element vertices */
        for (s = 0; s < nvert; s++)
            elem[s] = (int)elements[e+(s+1)*nelements];
			
		if (nvert == 3)
		{
			for (j = 0; j < 3; j++)
				for (s = 0; s < 3; s++)
					nod[s+3*j] = nodes[elem[s]+j*nnodes];
			inverse_matrix_tri(nod, accelerators[e].gradN);
		}
		
			
        for (j = 0; j < 3; j++)
        {
            accelerators[e].center[j] = 0.0;
            for (s = 0; s < nvert; s++)
                accelerators[e].center[j] += nodes[elem[s]+j*nnodes];
            accelerators[e].center[j] /= (double)nvert;
        }

        if (nvert == 3)
        {
            for (j = 0; j < 3; j++)
            {
                a[j] = nodes[elem[1]+j*nnodes]-nodes[elem[0]+j*nnodes];
                b[j] = nodes[elem[2]+j*nnodes]-nodes[elem[0]+j*nnodes];
            }
            cross(a, b, accelerators[e].n0);
        }
        else
        {
            for (j = 0; j < 3; j++)
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
                       const double *nodes,
                       int nelements,
                       const double *elements,
                       accelerator2D_t *accelerators)
{
#define NDIM 2
#define NVERT 2
    int e, j, s, nvert;
    int elem[NVERT];
	double a[NDIM];

    /* for each element */
    for (e = 0; e < nelements; e++)
    {
        for (s = 0; s < NVERT; s++)
            elem[s] = (int)elements[e+s*nelements];
			
        for (j = 0; j < NDIM; j++)
        {
            accelerators[e].center[j] = 0.0;
            for (s = 0; s < NVERT; s++)
                accelerators[e].center[j] += nodes[elem[s]+j*nnodes];
            accelerators[e].center[j] /= (double)nvert;
        }
		
		for (j = 0; j < NDIM; j++)
			a[j] = nodes[elem[1]+j*nnodes]-nodes[elem[0]+j*nnodes];

		accelerators[e].n0[0] = a[1];
		accelerators[e].n0[1] = -a[0];
    }
#undef NVERT
#undef NDIM
}

/* ------------------------------------------------------------------------ */
/* Determine Gaussian integration density based on distance between source */
/* and element center */
int gauss_division(const double *q,
                   const double *elemcenter,
                   const double *dist)
{
#define NDIM 3
    int j;
    double distance[NDIM], d;

    for (j = 0; j < NDIM; j++)
        distance[j] = elemcenter[j] - q[j];
    d = sqrt(dot(distance,distance));
    if (d < dist[0])
        return 1;
    if (d < dist[1])
        return 2;
    return 3;
#undef NDIM
}

/* ------------------------------------------------------------------------ */
/* Determine Gaussian integration density based on distance between source */
/* and element center */
int gauss_division2D(const double *q,
                   const double *elemcenter,
                   const double *dist)
{
#define NDIM 2
    int j;
    double distance[NDIM], d;

    for (j = 0; j < NDIM; j++)
        distance[j] = elemcenter[j] - q[j];
    d = sqrt(dot(distance,distance));
    if (d < dist[0])
        return 1;
    if (d < dist[1])
        return 2;
    return 3;
#undef NDIM
}
