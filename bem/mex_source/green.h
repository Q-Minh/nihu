#ifndef GREEN_H
#define GREEN_H

void green(const double *r,
           double k,
           double *gr,
           double *gi,
           const double *n,
           double *dgr,
           double *dgi);
		   
void green2D(const double *r,
           double k,
           double *gr,
           double *gi,
           const double *n,
           double *dgr,
           double *dgi);
		   
/* Second derivative of green function */
void ddgreen(const double *r,
	double k,
	const double *nx,
	const double *ny,
	double *ddgr,
	double *ddgi);

#endif
