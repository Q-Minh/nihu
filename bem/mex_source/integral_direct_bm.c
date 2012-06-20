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
/* Regular integral over a linear QUAD element                              */
void int_quad_lin_bm(const gauss_t *g,
                  const double *nodes,
                  const accelerator_t *accelerator,
                  const double *q, 						/* source point */
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
void int_quad_lin_sing_bm(const gauss_t *g,
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