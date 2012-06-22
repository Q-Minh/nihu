/* Burton-Miller calculation */
/* Intends to work on linear quad elements at the moment */
/* Peter Rucz, 2012-06-20 */

#include "integral_direct_bm.h"

#include "green.h"
#include "element.h"
#include "mesh.h"
#include "quadrature.h"

#include "vector.h" /* dot, cross */

#define _USE_MATH_DEFINES
#include <math.h>   /* sqrt, M_PI */
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
    int i, s;
	
    for (s = 0; s < 4; s++)
        ar[s] = ai[s] = br[s] = bi[s] = 0.0;

    /* for each gaussian integration point */
    for (i = 0; i < g->num; i++)
    {
		int j;
		double r[3], rxi[3], reta[3], norm[3], jac;
		double gr, gi, dgxr, dgxi, dgyr, dgyi, ddgr, ddgi;
		double Hr, Hi, Gr, Gi; /*Matrix elements */
		
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

		ddgreen(r, k, nq, norm, &gr, &gi, &dgxr, &dgxi, &dgyr, &dgyi, &ddgr, &ddgi);

		Hr = dgyr + alphar * ddgr - alphai * ddgi;
		Hi = dgyi + alphar * ddgi + alphai * ddgr;
		
		Gr = gr + alphar * dgxr - alphai * dgxi;
		Gi = gi + alphar * dgxi + alphai * dgxr;
		
		Hr *= jac;
		Hi *= jac;
		Gr *= jac;
		Gi *= jac;
		
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
/* Singular integral over a linear QUAD element - Burton-Miller             */
void int_quad_lin_sing_bm(const gauss_t *g, 	/* gaussian integration points and weights */
                       const double *nodes, 	/* element corner nodes */
                       const accelerator_t *accelerator,
                       int corner, 				/* singular corner index */
					   const double *nq, 		/* Normal vector at source point */
                       double k, 				/* wave number */
					   double alphar, 			/* real part of coupling constant */
					   double alphai,
                       double *ar, 				/* real part of matrix H or A */
                       double *ai,
                       double *br, 				/* real part of matrix  G  or B */
                       double *bi)
{
	double gradNq[4*3], Nq[4];
	/* Evaluate nabla N(x) and its normal component */
	{
		int s;
		double xiq, etaq;
		
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
		{
			double gradNqdotnx = dot(nq, gradNq+3*s) / 2.0;
			ar[s] = -alphar * gradNqdotnx;
			ai[s] = -alphai * gradNqdotnx;
			br[s] = bi[s] = 0.0;
			/* ar[s] = ai[s] = 0; */
		}
	}
	
	
	{
		int i, gnum;

		/* Preallocate for the singular quadrature */
		double *xiprime = (double *)malloc(sizeof(double)*2*g->num);
		double *etaprime = (double *)malloc(sizeof(double)*2*g->num);
		double *wprime = (double *)malloc(sizeof(double)*2*g->num);

		/* Obtain the singular quadrature */
		sing_quadr_corner_quad(g, corner, &gnum, xiprime, etaprime, wprime);

		/* for each gaussian integration point */
		for (i = 0; i < gnum; i++)
		{
			double N[4];
			int s, j;
			double r[3], norm[3], jac;
			
			double gr, gi, g0;
			double dgxr, dgxi, dgyr, dgyi, dg0x, dg0y;
			double ddgr, ddgi, ddg0;
			
			double gradNqdotny[4];
			
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

			/* Evaluate Green's functions */
			ddgreen(r, k, nq, norm, &gr, &gi, &dgxr, &dgxi, &dgyr, &dgyi, &ddgr, &ddgi);
			ddgreen0(r, nq, norm, &g0, &dg0x, &dg0y, &ddg0);
			
			for (s = 0; s < 4; ++s)
			{
				/* Term 1 */
				double tr = N[s] * k*k/8/M_PI/sqrt(dot(r,r));
				double ti = N[s] * k*k*k/12/M_PI;
				/* Term 2 */
				tr += (N[s] - Nq[s] - dot(gradNq+3*s, r)) * ddg0;
				/* Term 3 */
				tr += dot(gradNq+3*s, norm) * dg0x;
				
				/* H matrix: Add terms */
				ar[s] += (N[s] * dgyr + alphar * tr - alphai * ti) * jac;
				ai[s] += (N[s] * dgyi + alphar * ti + alphai * tr) * jac;
				
				br[s] += N[s]*(gr + alphar * dgxr - alphai * dgxi)  * jac;
				bi[s] += N[s]*(gi + alphar * dgxi + alphai * dgxr)  * jac;
			}
		}

		/* Free singular quadrature points */
		free(xiprime);
		free(etaprime);
		free(wprime);
	}
}
