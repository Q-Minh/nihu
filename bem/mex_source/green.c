#include "green.h"

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
    *gr = cos(k*ar) / ar;
    *gi = -sin(k*ar) / ar;
    if (n)
    {
        rn = dot(r, n) / ar;
        *dgr = (- *gr / ar + *gi * k) * rn;
        *dgi = (- *gi / ar - *gr * k) * rn;
    }
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
	
	gr = cos(k*ar)/ar;
	gi = -sin(k*ar)/ar;
	
	drnxdrny = -dot(r, nx)*dot(r, ny)/ar2;
	nxny = dot(nx, ny);
	
	br = (3.0/ar2 - k*k) * drnxdrny + nxny / ar2;
	bi = k/ar * (3.0*drnxdrny + nxny);
	
	*ddgr = gr*br - gi*bi;
	*ddgi = gr*bi + gi*br;
}
