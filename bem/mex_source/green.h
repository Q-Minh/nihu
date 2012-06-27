#ifndef GREEN_H
#define GREEN_H

/* exp(-1ikr) */
void greenr(const double *r,
			double k,
			double *gr,
			double *gi);

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
		   
/* 3D Green's function and its derivatives  */
void ddgreen(const double *r, /* y - x */
	double k,
	const double *nx,
	const double *ny,
	double *gr, 	/* Green's function */
	double *gi,
	double *dgxr, 	/* normal derivative with respect to x */
	double *dgxi,
	double *dgyr,	/* normal derivative with respect to y */
	double *dgyi,
	double *ddgr, 	/* double normal derivative */
	double *ddgi);

/* 3D static Green's function and its derivatives  */
void ddgreen0(const double *r, /* y - x */
	const double *nx,
	const double *ny,
	double *g, 	/* Green's function */
	double *dgx, 	/* normal derivative with respect to x */
	double *dgy,	/* normal derivative with respect to y */
	double *ddg); 	/* double normal derivative */
	
#endif
