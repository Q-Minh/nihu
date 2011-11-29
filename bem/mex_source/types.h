#ifndef TYPES_H
#define TYPES_H

// Structure containing Gaussian integration data
typedef struct
{
    int    num;		// number of Gaussian points
    double *N;		// Shape funtion values
    double *Nxi;	// Shape function derivatives
    double *Neta;
    double *w;		// Gaussian weights
    double *xi;		// Gaussian locations
} gauss_t;

#define SINGULAR_QUADRATURE_START 3

// Structure containing accelerator data for fast integration
typedef struct
{
	double center[3];	// element center
	double n0[3];		// normal at element center
	double nxi[3];		// derivatives of normal
	double neta[3];
	double gradN[3*3];
} accelerator_t;

typedef unsigned char boolean;

#endif
