#include "green.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include "vector.h"

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
    double ar, rn;

    ar = sqrt(dot(r, r));
    *gr = cos(k*ar) / ar / (4.0 * M_PI);
    *gi = -sin(k*ar) / ar / (4.0 * M_PI);
    if (n)
    {
        rn = dot(r, n) / ar;
        *dgr = (- *gr / ar + *gi * k) * rn;
        *dgi = (- *gi / ar - *gr * k) * rn;
    }
}

/* Compute Green's function and its second derivative (d^2 G / (d nx d ny)) */
void ddgreen(const double *r,		/* Location */
	double k, 						/* Wave number */
	const double* nx, 				/* n_x normal */
	const double* ny,				/* n_y normal */
	double *ddgr, 					/* real result */
	double *ddgi) 					/* imag result */
{
	double ar, ar2, ar3, rnx, rny, nxdny, ddr, r1, i1, r2, i2;
	ar2 = dot(r,r);
	ar = sqrt(ar2);								/* absolute distance */
	ar3 = ar2*ar;
    nxdny = dot(nx, ny); 						/* n_x dot n_y */
	ddr = dot(r, nx)*dot(r, ny) / ar2;			/* dr / (dn_x)* dr /dn_y */
	
	r1 = (3 - k*k*ar2)*ddr + nxdny; 			/* real part in brackets */
	i1 = (3*k*ar)*ddr + k*ar*nxdny; 			/* imaginary part in brackets */
	
	r2 =  cos(k*ar) / ar3 / (4.0 * M_PI);
	i2 = -sin(k*ar) / ar3 / (4.0 * M_PI);
	
	*ddgr = r1*r2 - i1*i2;
	*ddgi = r1*i2 + i1*r2;
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

void ddgreen(const double *r,
	double k,
	const double *nx,
	const double *ny,
	double *ddgr,
	double *ddgi)
{
	double ar, ar2, br, bi;
	double drnxdrny;
	double nxny;
	double gr, gi;
	
	ar2 = dot(r,r);
	ar = sqrt(ar2);
	
	gr = cos(k*ar)/ar / (4.0 * M_PI);
	gi = -sin(k*ar)/ar / (4.0 * M_PI);
	
	drnxdrny = -dot(r, nx)*dot(r, ny)/ar2;
	nxny = dot(nx, ny);
	
	br = (3.0/ar2 - k*k) * drnxdrny + nxny / ar2;
	bi = k/ar * (3.0*drnxdrny + nxny);
	
	*ddgr = gr*br - gi*bi;
	*ddgi = gr*bi + gi*br;
}
