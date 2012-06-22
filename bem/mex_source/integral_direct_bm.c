/* Burton-Miller calculation */
/* Intends to work on linear quad elements at the moment */
/* Peter Rucz, 2012-06-20 */

#include "integral_direct_bm.h"

#include "green.h"
#include "element.h"
#include "quadrature.h"

#include "vector.h" /* dot, cross */

#include <math.h>   /* sqrt */
#include <stdlib.h> /* malloc */


/* ------------------------------------------------------------------------ */
/* Regular integral over a linear QUAD element - Burton-Miller              */
void int_quad_lin_bm(const gauss_t *g,
					const double *nodes,
					const accelerator_t *accelerator,
					const double *q, 						/* source point */
					const double *nq, 						/* normal at source point */
					double k,
					double alphar,
					double alphai, 
                  double *ar,
                  double *ai,
                  double *br,
                  double *bi)
{
    int i, j, s;
    double r[3], rxi[3], reta[3], norm[3], jac, gr, gi, dgr, dgi, ddgr, ddgi;
	double Hr, Hi, Gr, Gi; /*Matrix elements */
	double grx, gix, dgrx, dgix; /* For x stuff */
	
    for (s = 0; s < 4; s++)
        ar[s] = ai[s] = br[s] = bi[s] = 0.0;

    /* for each gaussian integration point */
    for (i = 0; i < g->num; i++)
    {
        /* for each coordinate direction */
        for (j = 0; j < 3; j++)
        {
            r[j] = -q[j];
            rxi[j] = reta[j] = 0.0;
            /* computing integration location and its derivatives */
            for (s = 0; s < 4; s++)
            {
                r[j] += g->N[i+s*g->num]*nodes[s+4*j];
                rxi[j] += g->Nxi[i+s*g->num]*nodes[s+4*j];
                reta[j] += g->Neta[i+s*g->num]*nodes[s+4*j];
            }
        }
        /* surface normal and jacobian */
        cross(rxi, reta, norm);
        jac = sqrt(dot(norm, norm));
        for (j = 0; j < 3; j++)
            norm[j] /= jac;
        jac *= g->w[i];

        green(r, k, &gr, &gi, norm, &dgr, &dgi);
		green(r, k, &grx, &gix, nq, &dgrx, &dgix);
		ddgreen(r, k, norm, nq, &ddgr, &ddgi);

        gr = gr*jac;
        gi = gi*jac;
        dgr = dgr*jac;
        dgi = dgi*jac;
		dgrx = -dgrx*jac;
        dgix = -dgix*jac;
		ddgr = ddgr*jac;
		ddgi = ddgi*jac;

		Hr = dgr - alphar * ddgr + alphai * ddgi;
		Hi = dgi - alphar * ddgi - alphai * ddgr;
		
		Gr = gr - alphar * dgrx + alphai * dgix;
		Gi = gi - alphar * dgix - alphai * dgrx;
		
        for (s = 0; s < 4; s++)
        {
			/* H matrix */
            ar[s] += g->N[i+s*g->num]*Hr;
            ai[s] += g->N[i+s*g->num]*Hi;
			/* G matrix */
            br[s] += g->N[i+s*g->num]*Gr;
            bi[s] += g->N[i+s*g->num]*Gi;
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Singular integral over a linear QUAD element                             */
void int_quad_lin_sing_bm(const gauss_t *g,
                       const double *nodes,
                       const accelerator_t *accelerator,
                       int corner,
					   double *nq, 							/* Normal vector at source point */
                       double k,
					   double alphar,
					   double alphai,
                       double *ar,
                       double *ai,
                       double *br,
                       double *bi)
{
	
    int i, j, s, gnum;
    double r[3], norm[3], jac, gr, gi, dgr, dgi;
	
	double g0rx, g0ix;
	double ddgr, ddgi, ddg0r, ddg0i;
	double dgrx, dgix, dg0rx, dg0ix;
	
	double Nq[4], gradNq[12];
		
	double xiq, etaq;
	
	double dNdotny[4], dNdotnx[4];
	
	/* Terms in integral */
	double t1r, t1i;		/* Term 1: int N_j(y) [ d^2 Gk - d^2G0] */
	double t2, t2r, t2i;  		/* Term 2: (gradN(x) dot ny) * (dG0/dnx) */
	
	double t4[4]; 				/* Term 4: -1/2 grad N(x) * nx */ 
	
	/* Evaluate the shape functions */
	switch(corner)
	{
		case 0: xiq = -1.0; etaq = -1.0; break;
		case 1: xiq =  1.0; etaq = -1.0; break;
		case 2: xiq =  1.0; etaq =  1.0; break;
		case 3: xiq = -1.0; etaq =  1.0; break;
	}
	
	/* Evaluate shape functions and gradient */
	shapefun_quad(xiq, etaq, Nq);
	inverse_matrix_quad(nodes, xiq, etaq, gradNq);
	
	/* Initialize dNdotnx */
	for (s = 0; s < 4; s++)
		dNdotnx[s] = dot(nq, gradNq+3*s);
	
    double *xiprime, *etaprime, *wprime;
    double N[4];

	/* Initialize the result */
    for (s = 0; s < 4; s++)
        ar[s] = ai[s] = br[s] = bi[s] = 0.0;
	
	/* Calculate term4 */
	for (s = 0; s < 4; s++)
		t4[s] = -0.5 * dNdotnx[s];

	/* Preallocate for the singular quadrature */
    xiprime = (double *)malloc(sizeof(double)*2*g->num);
    etaprime = (double *)malloc(sizeof(double)*2*g->num);
    wprime = (double *)malloc(sizeof(double)*2*g->num);

	/* Obtain the singular quadrature */
    sing_quadr_corner_quad(g, corner, &gnum, xiprime, etaprime, wprime);

	/* TODO: obtain normal at x node */
	/* NOTE: taken as parameter */
	
    /* for each gaussian integration point */
    for (i = 0; i < gnum; i++)
    {
        shapefun_quad(xiprime[i], etaprime[i], N);
        /* computing integration location */
        for (j = 0; j < 3; j++)
        {
            r[j] = -nodes[corner+4*j];
            for (s = 0; s < 4; s++)
                r[j] += N[s]*nodes[s+4*j];

            norm[j] = accelerator->n0[j] + accelerator->nxi[j] * xiprime[i] + accelerator->neta[j] * etaprime[i];
        }

        jac = sqrt(dot(norm, norm));
        for (j = 0; j < 3; j++)
            norm[j] /= jac;
        jac *= wprime[i];

		/* Evaluate green function */
        green(r, k, &gr, &gi, norm, &dgr, &dgi);
		
		/* Evaluate hypersingular parts */
		/* #1 component */
		ddgreen(r, k, norm, nq, &ddgr, &ddgi);
		ddgreen0(r, k, norm, nq, &ddg0r, &ddg0i);
		
		green(r, k, &gr, &gi, norm, &dgr, &dgi);
		/* NOTE: sign for switching source and receiver? */
		green(r, k, &grx, &gix, nq, &dgrx, &dgix);
		green0(r, k, &g0rx, &g0ix, nq, &dg0rx, &dg0ix);

        gr *= jac;
        gi *= jac;
        dgr *= jac;
        dgi *= jac;
		ddgr *= jac;
		ddgi *= jac;
		dgrx *= jac;
		dgix *= jac;
		ddg0r *= jac;
		ddg0i *= jac;
		
		/* Term 1 */
		t1r = -(ddgr - ddg0r)*alphar + (ddgi - ddg0i)*alphai; t1r *= jac;
		t1i = -(ddgi - ddg0i)*alphar - (ddgr - ddg0r)*alphai; t1i *= jac;

        for (s = 0; s < 4; s++)
        {
			/* Calculate grad N * ny */
			dNdotny[s] = dot(norm, gradNq+3*s);
			
			/* Calculate term2 */
			t2 = dNdotny[s]*dg0rx;
			t2i = alphai*t2;
			t2r = alphar*t2;
			
			/* H matrix: Add term 1 */
			ar[s] += N[s]*t1r;
            ai[s] += N[s]*t1i;
			/* H matrix: Add term 4 */
			ar[s] += t4[s]*jac;
            			
            ar[s] += N[s]*dgr;
            ai[s] += N[s]*dgi;
            br[s] += N[s]*gr;
            bi[s] += N[s]*gi;
        }
    }

	/* Free singular quadrature points */
    free(xiprime);
    free(etaprime);
    free(wprime);
	
}