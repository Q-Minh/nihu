#include "green.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include "vector.h"

/* ------------------------------------------------------------------------ */
/* Compute Green's function and its normal derivative                       */
void greenr(const double *r,
           double k,
           double *gr,
           double *gi)
{
    double ar = sqrt(dot(r, r));
	/* exp(-ikr)/(4*pi*k)*i */
    *gi = cos(k*ar) / (4.0 * M_PI * k);
    *gr = sin(k*ar) / (4.0 * M_PI * k);
}

/* ------------------------------------------------------------------------ */
/* Compute Green's function and its normal derivative                       */
void green(const double *r,
           double k,
           double *gr,
           double *gi,
           const double *n,
           double *dgr,
           double *dgi)
{
    double ar = sqrt(dot(r, r));
    *gr = cos(k*ar) / ar / (4.0 * M_PI);
    *gi = -sin(k*ar) / ar / (4.0 * M_PI);
    if (n)
    {
        double rn = dot(r, n) / ar;
        *dgr = (- *gr / ar + *gi * k) * rn;
        *dgi = (- *gi / ar - *gr * k) * rn;
    }
}

/* ------------------------------------------------------------------------ */
/* Compute Green's function and its normal derivative                       */
void green2D(const double *r,
           double k,
           double *gr,
           double *gi,
           const double *n,
           double *dgr,
           double *dgi)
{
}

/* 3D Green's function and its derivatives  */
void ddgreen(const double *r, /* y - x */
	double k, 			/* wave number */
	const double *nx, 	/* unit normal at x */
	const double *ny,	/* unit normal at y */
	double *gr, 		/* Green's function */
	double *gi,
	double *dgxr, 		/* normal derivative with respect to x */
	double *dgxi,
	double *dgyr,		/* normal derivative with respect to y */
	double *dgyi,
	double *ddgr, 		/* double normal derivative */
	double *ddgi)
{
	double ar2 = dot(r,r); 			/* square norm of distance vector */
	double ar = sqrt(ar2); 			/* norm of distance vector */
    double rnx = -dot(r, nx) / ar;	/* normal derivative of distance wrt x */
    double rny = dot(r, ny) / ar;	/* normal derivative of distance wrt y */
	double rnxrny = rnx*rny;
	
	/* Green's function */
	*gr = cos(k*ar)/ar / (4.0 * M_PI);
	*gi = -sin(k*ar)/ar / (4.0 * M_PI);
	
	/* Normal derivatives */
    *dgxr = (- *gr / ar + *gi * k) * rnx;
    *dgxi = (- *gi / ar - *gr * k) * rnx;
    *dgyr = (- *gr / ar + *gi * k) * rny;
    *dgyi = (- *gi / ar - *gr * k) * rny;

	{
		double nxny = dot(nx, ny);
		double br = (3.0/ar2 - k*k) * rnxrny + nxny / ar2;
		double bi = k/ar * (3.0*rnxrny + nxny);
		/* Double normal derivative */
		*ddgr = *gr * br - *gi * bi;
		*ddgi = *gr * bi + *gi * br;
	}
}

/* 3D static Green's function and its derivatives  */
void ddgreen0(const double *r, /* y - x */
	const double *nx, 	/* unit normal at x */
	const double *ny, 	/* unit normal at y */
	double *g, 	/* Green's function */
	double *dgx, 	/* normal derivative with respect to x */
	double *dgy,	/* normal derivative with respect to y */
	double *ddg) 	/* double normal derivative */
{
	double ar2 = dot(r,r);
	double ar = sqrt(ar2);
	/* Normal derivatives */
    double rnx = -dot(r, nx) / ar;
    double rny = dot(r, ny) / ar;
	
	/* Green's function */
	*g = 1.0/(4.0 * M_PI * ar);
	
    *dgx = -*g / ar * rnx;
    *dgy = -*g / ar * rny;
	
	/* Double normal derivative */
	*ddg = *g/ar2 * (3.0 * rnx*rny + dot(nx, ny));
}
