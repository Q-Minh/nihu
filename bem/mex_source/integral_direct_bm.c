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
	double xi[16], w[16];
	int gnum = 16;
	/*
	xi[0] =  -0.906179845938664; w[0] = 0.236926885056189;
	xi[1] =  -0.538469310105683; w[1] = 0.478628670499367;
	xi[2] =   0.000000000000000; w[2] = 0.568888888888889;
	xi[3] =   0.538469310105683; w[3] = 0.478628670499367;
	xi[4] =   0.906179845938664; w[4] = 0.236926885056189;
	*/
	/*
	xi[0] = -0.968160239507626; w[0] = 0.081274388361575;
    xi[1] = -0.836031107326636; w[1] = 0.180648160694858;
    xi[2] = -0.613371432700590; w[2] = 0.260610696402935;
    xi[3] = -0.324253423403809; w[3] = 0.312347077040002;
    xi[4] =  0.000000000000000; w[4] = 0.330239355001259;
    xi[5] =  0.324253423403809; w[5] = 0.312347077040002;
    xi[6] =  0.613371432700591; w[6] = 0.260610696402936;
    xi[7] =  0.836031107326636; w[7] = 0.180648160694857;
    xi[8] =  0.968160239507626; w[8] = 0.081274388361575;
	*/
	
xi[ 0] = -0.9894009349916495;   w[ 0] = 0.0271524594117546;
xi[ 1] = -0.9445750230732327;   w[ 1] = 0.0622535239386470;
xi[ 2] = -0.8656312023878319;   w[ 2] = 0.0951585116824926;
xi[ 3] = -0.7554044083550033;   w[ 3] = 0.1246289712555342;
xi[ 4] = -0.6178762444026443;   w[ 4] = 0.1495959888165764;
xi[ 5] = -0.4580167776572275;   w[ 5] = 0.1691565193950028;
xi[ 6] = -0.2816035507792588;   w[ 6] = 0.1826034150449246;
xi[ 7] = -0.0950125098376372;   w[ 7] = 0.1894506104550683;
xi[ 8] = 0.0950125098376372;   w[ 8] = 0.1894506104550689;
xi[ 9] = 0.2816035507792590;   w[ 9] = 0.1826034150449239;
xi[10] = 0.4580167776572271;   w[10] = 0.1691565193950028;
xi[11] = 0.6178762444026440;   w[11] = 0.1495959888165774;
xi[12] = 0.7554044083550031;   w[12] = 0.1246289712555348;
xi[13] = 0.8656312023878318;   w[13] = 0.0951585116824925;
xi[14] = 0.9445750230732326;   w[14] = 0.0622535239386481;
xi[15] = 0.9894009349916498;   w[15] = 0.0271524594117544;
	
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
