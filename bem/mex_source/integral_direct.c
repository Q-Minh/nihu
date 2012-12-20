#include "integral_direct.h"

#include "green.h"
#include "element.h"
#include "quadrature.h"

#include "vector.h" /* dot, cross */

#include <math.h>   /* sqrt */
#include <stdlib.h> /* malloc */

/* ------------------------------------------------------------------------ */
/* Regular integral over a linear LINE element                              */
void int_line_lin(const gauss2D_t *g,
        const double *nodes,
        const accelerator2D_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    enum {NDIM = 2, NVERT = 2};
    int i, j, s;
    double r[NDIM], norm[NDIM], jac, gr, gi, dgr, dgi;
    
    for (s = 0; s < NVERT; s++)
        ar[s] = ai[s] = br[s] = bi[s] = 0.0;
    
    jac = sqrt(dot2D(accelerator->n0, accelerator->n0));
    for (j = 0; j < NDIM; j++)
        norm[j] = accelerator->n0[j]/jac;
    
    /* for each gaussian integration point */
    for (i = 0; i < g->num; i++)
    {
        /* computing integration location */
        for (j = 0; j < NDIM; j++)
        {
            r[j] = -q[j];
            for (s = 0; s < NVERT; s++)
                r[j] += g->N[i+s*g->num]*nodes[s+NDIM*j];
        }
        
        green2D(r, k, &gr, &gi, norm, &dgr, &dgi);
        
        gr *= g->w[i];
        gi *= g->w[i];
        dgr *= g->w[i];
        dgi *= g->w[i];
        
        for (s = 0; s < NVERT; s++)
        {
            ar[s] += g->N[i+s*g->num]*dgr;
            ai[s] += g->N[i+s*g->num]*dgi;
            br[s] += g->N[i+s*g->num]*gr;
            bi[s] += g->N[i+s*g->num]*gi;
        }
    }
    
    for (s = 0; s < NVERT; s++)
    {
        ar[s] *= jac;
        ai[s] *= jac;
        br[s] *= jac;
        bi[s] *= jac;
    }
}

/* ------------------------------------------------------------------------ */
/* Singular integral over a linear LINE element                             */
void int_line_lin_sing(const gauss2D_t *g,
        const double *nodes,
        const accelerator2D_t *accelerator,
        int corner,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    enum {NDIM = 2, NVERT = 2};
    int s;
    
    for (s = 0; s < NVERT; s++)
        ar[s] = ai[s] = br[s] = bi[s] = 0.0;
}

/* ------------------------------------------------------------------------ */
/* Regular integral over a constant LINE element                            */
void int_line_const(const gauss2D_t *g,
        const double *nodes,
        const accelerator2D_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    enum {NDIM = 2, NVERT = 2};
    int i, j, s;
    double r[NDIM], norm[NDIM], jac, gr, gi, dgr, dgi;
    
    *ar = *ai = *br = *bi = 0.0;
    
    jac = sqrt(dot2D(accelerator->n0, accelerator->n0));
    for (j = 0; j < NDIM; j++)
        norm[j] = accelerator->n0[j]/jac;
    
    /* for each gaussian integration point */
    for (i = 0; i < g->num; i++)
    {
        /* computing integration location */
        for (j = 0; j < NDIM; j++)
        {
            r[j] = -q[j];
            for (s = 0; s < NVERT; s++)
                r[j] += g->N[i+s*g->num]*nodes[s+NDIM*j];
        }
        
        green2D(r, k, &gr, &gi, norm, &dgr, &dgi);
        
        *ar += dgr*g->w[i];
        *ai += dgi*g->w[i];
        *br += gr*g->w[i];
        *bi += gi*g->w[i];
    }
    
    *ar *= jac;
    *ai *= jac;
    *br *= jac;
    *bi *= jac;
}

/* ------------------------------------------------------------------------ */
/* Singular integral over a constant LINE element                           */
void int_line_const_sing(const gauss2D_t *g,
        const double *nodes,
        const accelerator2D_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    enum {NDIM = 2, NVERT = 2};
    *ar = *ai = *br = *bi = 0.0;
}

