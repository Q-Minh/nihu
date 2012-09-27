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
#define NDIM 2
#define NVERT 2
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
#undef NDIM
#undef NVERT
}

/* ------------------------------------------------------------------------ */
/* Singular integral over a linear TRIA element                             */
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
#define NDIM 2
#define NVERT 2
    int s;

for (s = 0; s < NVERT; s++)
    ar[s] = ai[s] = br[s] = bi[s] = 0.0;
#undef NDIM
#undef NVERT
}

/* ------------------------------------------------------------------------ */
/* Regular integral over a linear TRIA element                              */
void int_tri_lin(const gauss_t *g,
        const double *nodes,
        const accelerator_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
#define NDIM 3
#define NVERT 3
    int i, j, s;
double r[NDIM], norm[NDIM], jac, gr, gi, dgr, dgi;

for (s = 0; s < NVERT; s++)
    ar[s] = ai[s] = br[s] = bi[s] = 0.0;

jac = sqrt(dot(accelerator->n0, accelerator->n0));
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
    
    green(r, k, &gr, &gi, norm, &dgr, &dgi);
    
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
#undef NDIM
#undef NVERT
}

/* ------------------------------------------------------------------------ */
/* Singular integral over a linear TRIA element                             */
void int_tri_lin_sing(const gauss_t *g,
        const double *nodes,
        const accelerator_t *accelerator,
        int corner,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    int i, j, s, gnum;
    double r[3], norm[3], jac, gr, gi, dgr, dgi;
    double *xiprime, *etaprime, *wprime;
    double N[3];
    
    for (s = 0; s < 3; s++)
        ar[s] = ai[s] = br[s] = bi[s] = 0.0;
    
    jac = sqrt(dot(accelerator->n0, accelerator->n0));
    for (j = 0; j < 3; j++)
        norm[j] = accelerator->n0[j]/jac;
    
    xiprime = (double *)malloc(sizeof(double)*1*g->num);
    etaprime = (double *)malloc(sizeof(double)*1*g->num);
    wprime = (double *)malloc(sizeof(double)*1*g->num);
    
    sing_quadr_corner_tria(g, corner, &gnum, xiprime, etaprime, wprime);
    
    /* for each gaussian integration point */
    for (i = 0; i < gnum; i++)
    {
        shapefun_tria(xiprime[i], etaprime[i], N);
        
        /* computing integration location */
        for (j = 0; j < 3; j++)
        {
            r[j] = -nodes[corner+3*j];
            for (s = 0; s < 3; s++)
                r[j] += N[s]*nodes[s+3*j];
        }
        
        green(r, k, &gr, &gi, norm, &dgr, &dgi);
        
        gr *= wprime[i];
        gi *= wprime[i];
        dgr *= wprime[i];
        dgi *= wprime[i];
        
        for (s = 0; s < 3; s++)
        {
            ar[s] += N[s]*dgr;
            ai[s] += N[s]*dgi;
            br[s] += N[s]*gr;
            bi[s] += N[s]*gi;
        }
    }
    
    for (s = 0; s < 3; s++)
    {
        ar[s] *= jac;
        ai[s] *= jac;
        br[s] *= jac;
        bi[s] *= jac;
    }
    
    free(xiprime);
    free(etaprime);
    free(wprime);
}

/* ------------------------------------------------------------------------ */
/* Regular integral over a linear QUAD element                              */
void int_quad_lin(const gauss_t *g,
        const double *nodes,
        const accelerator_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    int i, j, s;
    double r[3], rxi[3], reta[3], norm[3], jac, gr, gi, dgr, dgi;
    
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
        
        gr = gr*jac;
        gi = gi*jac;
        dgr = dgr*jac;
        dgi = dgi*jac;
        
        for (s = 0; s < 4; s++)
        {
            ar[s] += g->N[i+s*g->num]*dgr;
            ai[s] += g->N[i+s*g->num]*dgi;
            br[s] += g->N[i+s*g->num]*gr;
            bi[s] += g->N[i+s*g->num]*gi;
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Singular integral over a linear QUAD element                             */
void int_quad_lin_sing(const gauss_t *g,
        const double *nodes,
        const accelerator_t *accelerator,
        int corner,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    int i, j, s, gnum;
    double r[3], norm[3], jac, gr, gi, dgr, dgi;
    double *xiprime, *etaprime, *wprime;
    double N[4];
    
    for (s = 0; s < 4; s++)
        ar[s] = ai[s] = br[s] = bi[s] = 0.0;
    
    xiprime = (double *)malloc(sizeof(double)*2*g->num);
    etaprime = (double *)malloc(sizeof(double)*2*g->num);
    wprime = (double *)malloc(sizeof(double)*2*g->num);
    
    sing_quadr_corner_quad(g, corner, &gnum, xiprime, etaprime, wprime);
    
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
        
        green(r, k, &gr, &gi, norm, &dgr, &dgi);
        
        gr *= jac;
        gi *= jac;
        dgr *= jac;
        dgi *= jac;
        
        for (s = 0; s < 4; s++)
        {
            ar[s] += N[s]*dgr;
            ai[s] += N[s]*dgi;
            br[s] += N[s]*gr;
            bi[s] += N[s]*gi;
        }
    }
    
    free(xiprime);
    free(etaprime);
    free(wprime);
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
#define NDIM 2
#define NVERT 2
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
#undef NDIM
#undef NVERT
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
#define NDIM 2
#define NVERT 2
    *ar = *ai = *br = *bi = 0.0;
#undef NDIM
#undef NVERT
}

/* ------------------------------------------------------------------------ */
/* Regular integral over a constant TRIA element                            */
void int_tri_const(const gauss_t *g,
        const double *nodes,
        const accelerator_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    int i, j, s;
    double r[3], norm[3], jac, gr, gi, dgr, dgi;
    
    *ar = *ai = *br = *bi = 0.0;
    
    jac = sqrt(dot(accelerator->n0, accelerator->n0));
    for (j = 0; j < 3; j++)
        norm[j] = accelerator->n0[j]/jac;
    
    /* for each gaussian integration point */
    for (i = 0; i < g->num; i++)
    {
        /* computing integration location */
        for (j = 0; j < 3; j++)
        {
            r[j] = -q[j];
            for (s = 0; s < 3; s++)
                r[j] += g->N[i+s*g->num]*nodes[s+3*j];
        }
        
        green(r, k, &gr, &gi, norm, &dgr, &dgi);
        
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
/* Singular integral over a constant TRIA element                           */
void int_tri_const_sing(const gauss_t *g,
        const double *nodes,
        const accelerator_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    int i, j, s;
    double r[3], norm[3], jac, gr, gi, dgr, dgi;
    double *xiprime, *etaprime, *wprime;
    double N[3];
    int gnum;
    
    *ar = *ai = *br = *bi = 0.0;
    
    jac = sqrt(dot(accelerator->n0, accelerator->n0));
    for (j = 0; j < 3; j++)
        norm[j] = accelerator->n0[j] / jac;
    
    xiprime = (double *)malloc(sizeof(double)*3*g->num);
    etaprime = (double *)malloc(sizeof(double)*3*g->num);
    wprime = (double *)malloc(sizeof(double)*3*g->num);
    
    sing_quadr_face_tria(g, 1.0/3.0, 1.0/3.0, &gnum, xiprime, etaprime, wprime);
    
    /* for each gaussian integration point */
    for (i = 0; i < gnum; i++)
    {
        shapefun_tria(xiprime[i], etaprime[i], N);
        /* for each coordinate direction */
        for (j = 0; j < 3; j++)
        {
            r[j] = -q[j];
            /* computing integration location and its derivatives */
            for (s = 0; s < 3; s++)
                r[j] += N[s]*nodes[s+3*j];
        }
        
        green(r, k, &gr, &gi, norm, &dgr, &dgi);
        
        *ar += dgr*wprime[i];
        *ai += dgi*wprime[i];
        *br += gr*wprime[i];
        *bi += gi*wprime[i];
    }
    
    *ar *= jac;
    *ai *= jac;
    *br *= jac;
    *bi *= jac;
    
    free(xiprime);
    free(etaprime);
    free(wprime);
}


/* ------------------------------------------------------------------------ */
/* Regular integral over a constant QUAD element                            */
void int_quad_const(const gauss_t *g,
        const double *nodes,
        const accelerator_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    int i;
    
    *ar = *ai = *br = *bi = 0.0;
    
    /* for each gaussian integration point */
    for (i = 0; i < g->num; i++)
    {
        int j, s;
        double r[3], rxi[3], reta[3], norm[3], jac, gr, gi, dgr, dgi;
        
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
        
        *ar += dgr*jac;
        *ai += dgi*jac;
        *br += gr*jac;
        *bi += gi*jac;
    }
}

/* ------------------------------------------------------------------------ */
/* Singular integral over a constant QUAD element                           */
void int_quad_const_sing(const gauss_t *g,
        const double *nodes,
        const accelerator_t *accelerator,
        const double *q,
        double k,
        double *ar,
        double *ai,
        double *br,
        double *bi)
{
    int i, j, s;
    double r[3], norm[3], jac, gr, gi, dgr, dgi;
    double *xiprime, *etaprime, *wprime;
    double N[4];
    int gnum;
    
    *ar = *ai = *br = *bi = 0.0;
    
    xiprime = (double *)malloc(sizeof(double)*4*g->num);
    etaprime = (double *)malloc(sizeof(double)*4*g->num);
    wprime = (double *)malloc(sizeof(double)*4*g->num);
    
    sing_quadr_face_quad(g, 0.0, 0.0, &gnum, xiprime, etaprime, wprime);
    
    /* for each gaussian integration point */
    for (i = 0; i < gnum; i++)
    {
        shapefun_quad(xiprime[i], etaprime[i], N);
        /* for each coordinate direction */
        for (j = 0; j < 3; j++)
        {
            norm[j] = accelerator->n0[j] + accelerator->nxi[j] * xiprime[i] + accelerator->neta[j]*etaprime[i];
            r[j] = -q[j];
            /* computing integration location and its derivatives */
            for (s = 0; s < 4; s++)
                r[j] += N[s]*nodes[s+4*j];
        }
        
        /* surface normal and jacobian */
        jac = sqrt(dot(norm, norm));
        for (j = 0; j < 3; j++)
            norm[j] /= jac;
        jac *= wprime[i];
        
        green(r, k, &gr, &gi, norm, &dgr, &dgi);
        
        *ar += dgr*jac;
        *ai += dgi*jac;
        *br += gr*jac;
        *bi += gi*jac;
    }
    
    free(xiprime);
    free(etaprime);
    free(wprime);
}
