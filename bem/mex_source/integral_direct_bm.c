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
#define NVERT 3
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
        for (s = 0; s < NVERT; s++) 			/* for all vertices */
            r[j] += g->N[i+s*g->num]*nodes[s+NVERT*j];
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
#undef NVERT
}


/* ------------------------------------------------------------------------ */
/* Regular integral over a constant TRIA element using Burton-Miller        */
void int_quad_const_bm(const gauss_t *g,
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
#define NVERT 4
    int i,j,s;

/* Initialize result to zero */
*ar = *ai = *br = *bi = 0.0;

/* for each gaussian integration point */
for (i = 0; i < g->num; i++)
{
    int j, s;
    double r[3], rxi[3], reta[3], norm[3], jac, gr, gi, dgxr, dgxi, dgyr, dgyi, ddgr, ddgi;
    
    /* for each coordinate direction */
    for (j = 0; j < 3; j++)
    {
        r[j] = -q[j];
        rxi[j] = reta[j] = 0.0;
        /* computing integration location and its derivatives */
        for (s = 0; s < NVERT; s++)
        {
            r[j] += g->N[i+s*g->num]*nodes[s+NVERT*j];
            rxi[j] += g->Nxi[i+s*g->num]*nodes[s+NVERT*j];
            reta[j] += g->Neta[i+s*g->num]*nodes[s+NVERT*j];
        }
    }
    
    /* surface normal and jacobian */
    cross(rxi, reta, norm);
    jac = sqrt(dot(norm, norm));
    for (j = 0; j < 3; j++)
        norm[j] /= jac;
    jac *= g->w[i];
    
    /* Evaluate Green function and its derivatives */
    ddgreen(r, k, nq, norm, &gr, &gi, &dgxr, &dgxi, &dgyr, &dgyi, &ddgr, &ddgi);
    
    /* Evaluate matrix elements */
    *ar += (dgyr + alphar * ddgr - alphai * ddgi)*jac;
    *ai += (dgyi + alphar * ddgi + alphai * ddgr)*jac;
    
    *br += (gr + alphar * dgxr - alphai * dgxi)*jac;
    *bi += (gi + alphar * dgxi + alphai * dgxr)*jac;
}

#undef NVERT
}


double gauss_xi_bm_sing[] = {
    -0.995187219997022,
    -0.974728555971310,
    -0.938274552002733,
    -0.886415527004401,
    -0.820001985973903,
    -0.740124191578554,
    -0.648093651936975,
    -0.545421471388839,
    -0.433793507626045,
    -0.315042679696164,
    -0.191118867473616,
    -0.064056892862606,
    0.064056892862605,
    0.191118867473617,
    0.315042679696164,
    0.433793507626045,
    0.545421471388839,
    0.648093651936975,
    0.740124191578555,
    0.820001985973903,
    0.886415527004401,
    0.938274552002733,
    0.974728555971310,
    0.995187219997021
};
double gauss_w_bm_sing[] = {
    0.012341229799987,
    0.028531388628934,
    0.044277438817420,
    0.059298584915436,
    0.073346481411081,
    0.086190161531953,
    0.097618652104114,
    0.107444270115965,
    0.115505668053726,
    0.121670472927803,
    0.125837456346828,
    0.127938195346752,
    0.127938195346752,
    0.125837456346828,
    0.121670472927803,
    0.115505668053726,
    0.107444270115965,
    0.097618652104114,
    0.086190161531953,
    0.073346481411080,
    0.059298584915437,
    0.044277438817420,
    0.028531388628933,
    0.012341229799988
};

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
#define NVERT 3
    int i, j, s;

/* TODO: temporary */
int gnum = sizeof(gauss_w_bm_sing)/sizeof(gauss_w_bm_sing[0]);

/* Initialize the result */
*ar = *ai = *br = *bi = 0.0;

/* Go through all three sides */
for (i = 0; i < NVERT; ++i)
{
    double L; 				     	/* Length of side */
    int n1, n2;						/* The two node numbers */
    int ig;
    double d[3];
    /* Initialize nodes */
    n1 = i; n2 = (n1+1)%NVERT; 			/* Calculate between ith and i+1th nodes */
    /* Obtain distance vector */
    for (j = 0; j < 3; ++j) /* For all dimensions */
        d[j] = nodes[n2+NVERT*j] - nodes[n1+NVERT*j];
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
            r[j] = 0.5*(1.0-gauss_xi_bm_sing[ig])*nodes[n1+NVERT*j]
                    + 0.5*(1.0+gauss_xi_bm_sing[ig])*nodes[n2+NVERT*j] - q[j];
        /* Absolute value of distance */
        lr = sqrt(dot(r,r));
        /* Jacobian is sin(beta)/ar*L/2 */
        /* Weight is also part of jacobian */
        tmp = dot(r,d)/(lr*L);
        jac = gauss_w_bm_sing[ig] * sqrt(1.0 - tmp*tmp) / lr * L / 2.0;
        
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
#undef NVERT
}


/* ------------------------------------------------------------------------ */
/* Singular integral over a constant QUAD element using Burton-Miller       */
void int_quad_const_sing_bm(const gauss_t *g,   /* This will use line gauss! */
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
#define NVERT 4
    int i, j, s;

/* TODO: temporary */
int gnum = sizeof(gauss_w_bm_sing)/sizeof(gauss_w_bm_sing[0]);

/* Initialize the result */
*ar = *ai = *br = *bi = 0.0;

/* Go through all three sides */
for (i = 0; i < NVERT; ++i)
{
    double L; 				     	/* Length of side */
    int n1, n2;						/* The two node numbers */
    int ig;
    double d[3];
    /* Initialize nodes */
    n1 = i; n2 = (n1+1)%NVERT; 			/* Calculate between ith and i+1th nodes */
    /* Obtain distance vector */
    for (j = 0; j < 3; ++j) /* For all dimensions */
        d[j] = nodes[n2+NVERT*j] - nodes[n1+NVERT*j];
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
            r[j] = 0.5*(1.0-gauss_xi_bm_sing[ig])*nodes[n1+NVERT*j]
                    + 0.5*(1.0+gauss_xi_bm_sing[ig])*nodes[n2+NVERT*j] - q[j];
        /* Absolute value of distance */
        lr = sqrt(dot(r,r));
        /* Jacobian is sin(beta)/ar*L/2 */
        /* Weight is also part of jacobian */
        tmp = dot(r,d)/(lr*L);
        jac = gauss_w_bm_sing[ig] * sqrt(1.0 - tmp*tmp) / lr * L / 2.0;
        
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
#undef NVERT
}
