/**
* \file types.h
* \brief basic type definitions
* \author Peter Fiala fiala@hit.bme.hu
*/

#ifndef TYPES_H
#define TYPES_H

/**
* \brief Gaussian integration data
*/
struct gauss_t
{
    int    num;		/**< number of Gaussian points */
    double *N;		/**< Shape funtion values */
    double *Nxi;	/**< Shape function derivatives */
    double *Neta;	/**< Shape function derivatives */
    double *w;		/**< Gaussian weights */
    double *xi;		/**< Gaussian locations */
};

/**
* \brief  Gaussian integration data for 2D integration
*/
struct gauss2D_t
{
    int    num;		/**< number of Gaussian points */
    double *N;		/**< Shape funtion values */
    double *Nxi;	/**< Shape function derivatives */
    double *w;		/**< Gaussian weights */
    double *xi;		/**< Gaussian locations */
};

/**
* \brief  Accelerator structure for fast integration
*/
struct accelerator_t
{
	double center[3];	/**< element center */
	double n0[3];		/**< normal at element center */
	double nxi[3];		/**< derivatives of normal */
	double neta[3];		/**< derivative of normal */
	double gradN[3*3];	/**< shape function gradient matrix */
};

/**
* \brief  Accelerator structure for fast integration
*/
struct accelerator2D_t
{
	double center[2];	/**< element center */
	double n0[2];		/**< element normal */
};

#endif

