#include "green.hpp"

#include "vector.h"

const complex_scalar compJ = complex_scalar(0.0, 1.0);

complex_scalar greenr(const double ar,
           const complex_scalar &k)
{
	/* exp(-ikr)/(4*pi*k) */
    return exp(-compJ * k*ar) / (4.0 * M_PI * k);
}


complex_scalar green(const double *r,
           const complex_scalar &k)
{
    double ar = sqrt(dot(r, r));
    return exp(-compJ*k*ar) / (4.0 * M_PI * ar);
}


void green(const double *r,
           const complex_scalar &k,
           complex_scalar &g,
           const double *n,
           complex_scalar &dg)
{
    double ar = sqrt(dot(r, r));
    g = exp(-compJ*k*ar) / (4.0 * M_PI * ar);
    double rn = dot(r, n) / ar;
    dg = -g * (1.0 + compJ*k*ar)/ar * rn;
}

void green2D(const double *r,
             const complex_scalar &k,
             complex_scalar &g,
             const double *n,
             complex_scalar &dg)
{
}

void ddgreen(const double *r, /* y - x */
             const complex_scalar &k,
             const double *nx, 	/* unit normal at x */
             const double *ny,	/* unit normal at y */
             complex_scalar &g, 		/* Green's function */
             complex_scalar &dgx, 		/* normal derivative with respect to x */
             complex_scalar &dgy,		/* normal derivative with respect to y */
             complex_scalar &ddg)
{
    /* Green's function */
    double ar2 = dot(r,r); 			/* square norm of distance vector */
    double ar = sqrt(ar2); 			/* norm of distance vector */
    g = exp(-compJ*k*ar) / (4.0 * M_PI * ar);

    /* Normal derivatives */
    double rnx = -dot(r, nx) / ar;	/* normal derivative of distance wrt x */
    double rny = dot(r, ny) / ar;		/* normal derivative of distance wrt y */

    dgx = -g * (1.0 + compJ*k*ar)/ar * rnx;
    dgy = -g * (1.0 + compJ*k*ar)/ar * rny;

    /* Double normal derivative */
    ddg = g * ((3.0/ar2 + 3.0*compJ*k/ar-k*k)*rnx*rny + 1.0/ar2*(1.0+compJ*k*ar)*dot(nx, ny));
}

/* 3D static Green's function and its derivatives  */
void ddgreen0(const double *r, /* y - x */
              const double *nx, 	/* unit normal at x */
              const double *ny, 	/* unit normal at y */
              double &g, 	/* Green's function */
              double &dgx, 	/* normal derivative with respect to x */
              double &dgy,	/* normal derivative with respect to y */
              double &ddg) 	/* double normal derivative */
{

    /* Green's function */
    double ar2 = dot(r,r);
    double ar = sqrt(ar2);
    g = 1.0/(4.0 * M_PI * ar);

    /* normal derivatives */
    double rnx = -dot(r, nx) / ar;
    double rny = dot(r, ny) / ar;
    dgx = -g / ar * rnx;
    dgy = -g / ar * rny;

    /* Double normal derivative */
    ddg = g/ar2 * (3.0 * rnx*rny + dot(nx,ny));
}
