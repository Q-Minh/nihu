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

/* ------------------------------------------------------------------------ */
/* Regular integral over a constant TRIA element using Burton-Miller        */
void int_tri_const_bm(const gauss_t *g,
					const double *nodes,
					const accelerator_t *accelerator,
					const double *q,					/* source location */
					const double *nq, 					/* source normal vector */
					double k, 							/* Wave number */
					double alphar,						/* Coupling real */
					double alphai,						/* Coupling imag */
					double *ar,
					double *ai,
					double *br,
					double *bi)
{
    int i,j,s;
    double norm[3], jac;

	/* Initialize result to zero */
    *ar = *ai = *br = *bi = 0.0;

	/* Jacobian and surface normal calculation  */
    jac = sqrt(dot(accelerator->n0, accelerator->n0));
    for (j = 0; j < 3; j++)
        norm[j] = accelerator->n0[j]/jac;

    /* for each gaussian integration point */
    for (i = 0; i < g->num; i++)
    {
		double r[3];
		double gr, gi, dgxr, dgxi, dgyr, dgyi, ddgr, ddgi;
		
        /* computing integration location */
        for (j = 0; j < 3; j++) 				/* for all directions */
        {
            r[j] = -q[j];
            for (s = 0; s < 3; s++) 			/* for all vertices */
                r[j] += g->N[i+s*g->num]*nodes[s+3*j];
        }

        /* Evaluate Green function and its derivatives */
        ddgreen(r, k, nq, norm, &gr, &gi, &dgxr, &dgxi, &dgyr, &dgyi, &ddgr, &ddgi);

		/* Evaluate matrix elements */
		*ar += (dgyr + alphar * ddgr - alphai * ddgi)*(g->w[i]);
		*ai += (dgyi + alphar * ddgi + alphai * ddgr)*(g->w[i]);
		
		*br += (gr + alphar * dgxr - alphai * dgxi)*(g->w[i]);
		*bi += (gi + alphar * dgxi + alphai * dgxr)*(g->w[i]);
    }

    /* Finally, multiply with jacobian */
    *ar *= jac;
    *ai *= jac;
    *br *= jac;
    *bi *= jac;
}

/* ------------------------------------------------------------------------ */
/* Singular integral over a constant TRIA element using Burton-Miller       */
void int_tri_const_sing_bm(const gauss_t *g,   /* This will use line gauss! */
						const double *nodes,
						const accelerator_t *accelerator,
						const double *q, 		/* Source location */
						const double *nq, 		/* Source normal */
						double k, 				/* Wave number */
						double alphar,			/* Coupling constant real */
						double alphai, 			/* Coupling constant imag */
						double *ar,
						double *ai,
						double *br,
						double *bi)
{
    int i, j, s;
	
	/* TODO: temporary */
    double xi[24], w[24];
	int gnum = 24;
	
xi[ 0] = -0.9951872199970215;   w[ 0] = 0.0123412297999869;
xi[ 1] = -0.9747285559713099;   w[ 1] = 0.0285313886289337;
xi[ 2] = -0.9382745520027330;   w[ 2] = 0.0442774388174201;
xi[ 3] = -0.8864155270044012;   w[ 3] = 0.0592985849154364;
xi[ 4] = -0.8200019859739028;   w[ 4] = 0.0733464814110806;
xi[ 5] = -0.7401241915785545;   w[ 5] = 0.0861901615319534;
xi[ 6] = -0.6480936519369753;   w[ 6] = 0.0976186521041143;
xi[ 7] = -0.5454214713888395;   w[ 7] = 0.1074442701159653;
xi[ 8] = -0.4337935076260454;   w[ 8] = 0.1155056680537255;
xi[ 9] = -0.3150426796961637;   w[ 9] = 0.1216704729278035;
xi[10] = -0.1911188674736164;   w[10] = 0.1258374563468276;
xi[11] = -0.0640568928626056;   w[11] = 0.1279381953467523;
xi[12] = 0.0640568928626054;   w[12] = 0.1279381953467524;
xi[13] = 0.1911188674736169;   w[13] = 0.1258374563468279;
xi[14] = 0.3150426796961636;   w[14] = 0.1216704729278030;
xi[15] = 0.4337935076260452;   w[15] = 0.1155056680537262;
xi[16] = 0.5454214713888393;   w[16] = 0.1074442701159654;
xi[17] = 0.6480936519369753;   w[17] = 0.0976186521041142;
xi[18] = 0.7401241915785546;   w[18] = 0.0861901615319534;
xi[19] = 0.8200019859739028;   w[19] = 0.0733464814110800;
xi[20] = 0.8864155270044011;   w[20] = 0.0592985849154368;
xi[21] = 0.9382745520027326;   w[21] = 0.0442774388174205;
xi[22] = 0.9747285559713097;   w[22] = 0.0285313886289333;
xi[23] = 0.9951872199970212;   w[23] = 0.0123412297999878;

	/* Initialize the result */
    *ar = *ai = *br = *bi = 0.0;
	
	/* Go through all three sides */
    for (i = 0; i < 3; ++i)
	{
		double L; 				     	/* Length of side */
		int n1, n2;						/* The two node numbers */
		int ig;
		double d[3];					
		/* Initialize nodes */
		n1 = i; n2 = (n1+1)%3; 			/* Calculate between ith and i+1th nodes */
		/* Obtain distance vector */
		for (j = 0; j < 3; ++j) /* For all dimensions */
			d[j] = nodes[n2+3*j] - nodes[n1+3*j];
		/* Calculate element length */
		L = sqrt(dot(d,d));
		
		/* Go through all integration points */
		for (ig = 0; ig < gnum; ig++)
		{
			double r[3], lr, jac;  		/* actual r vector, r, and jacobian */
			double gr, gi, grr, gri;
			double tmp;
			/* Calculate actual location x(\xi)-x_q*/
			for (j=0; j < 3; ++j)
				r[j] = 0.5*(1.0-xi[ig])*nodes[n1+3*j] + 0.5*(1.0+xi[ig])*nodes[n2+3*j] - q[j];
			/* Absolute value of distance */
			lr = sqrt(dot(r,r));
			/* Jacobian is sin(beta)/ar*L/2 */
			/* Weight is also part of jacobian */
			tmp = dot(r,d)/(lr*L);
			jac = w[ig] * sqrt(1.0 - tmp*tmp) / lr * L / 2.0;
			
			if (jac < 0)
				mexPrintf("Jac < 0!\n");
			
			/* Here the calculation of the integrand should be performed */
			/* Matrix H: the negative of the simple green function should be evaluated */
			green(r, k, &gr, &gi, NULL, NULL, NULL);

			*ar -= (gr * alphar - gi*alphai)*jac;
			*ai -= (gr * alphai + gi*alphar)*jac;
			
			/* Matrix G: the reduced Green is evaluated */
			greenr(lr, k, &grr, &gri);
			
			*br += -gri*jac;
			*bi += grr*jac;
			
		}
	}
	
}
