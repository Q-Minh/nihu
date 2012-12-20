#include "integral_indirect.h"

#include "green.h"
#include "element.h"
#include "quadrature.h"

#include "vector.h" /* dot, cross */

#include <math.h>   /* sqrt */
#include <stdlib.h> /* malloc */


/* ------------------------------------------------------------------------ */
/* Regular double integral over two constant TRIA elements                  */
void dint_tritri_const(const gauss_t *g3,
                       const double *nodesx,
                       const accelerator_t *acceleratorx,
                       const double *nodesy,
                       const accelerator_t *acceleratory,
                       double k,
                       double *br,
                       double *bi)
{
    int ix, iy, j, s;
    double rx[3], jacx, jacy;
    double r[3], gr, gi;

    *br = *bi = 0.0;

    jacx = sqrt(dot(acceleratorx->n0, acceleratorx->n0));
    jacy = sqrt(dot(acceleratory->n0, acceleratory->n0));

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each y gaussian integration point */
        for (iy = 0; iy < g3->num; iy++)
        {
            /* computing integration y location */
            for (j = 0; j < 3; j++)
            {
                r[j] = -rx[j];
                for (s = 0; s < 3; s++)
                    r[j] += g3->N[iy+s*g3->num]*nodesy[s+3*j];
            }

            green(r, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr*g3->w[ix]*g3->w[iy];
            *bi += gi*g3->w[ix]*g3->w[iy];
        }
    }

    *br *= jacx*jacy;
    *bi *= jacx*jacy;
}

/* ------------------------------------------------------------------------ */
/* Compute regular double B integral over two constant QUAD elements        */
void dint_quadquad_const(const gauss_t *g4,
                         const double *nodesx,
                         const accelerator_t *acceleratorx,
                         const double *nodesy,
                         const accelerator_t *acceleratory,
                         double k,
                         double *br,
                         double *bi)
{
    int ix, iy, j, s;
    double rx[3], ry[3], norm[3], jacx, jacy, xi, eta;
    double gr, gi;

    *br = *bi = 0.0;

    /* for each x gaussian integration point */
    for (ix = 0; ix < g4->num; ix++)
    {
        xi = g4->xi[ix];
        eta = g4->xi[ix+g4->num];
        /* for each coordinate direction */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 4; s++)
                rx[j] += g4->N[ix+s*g4->num]*nodesx[s+4*j];
            norm[j] = acceleratorx->n0[j] + xi*acceleratorx->nxi[j] + eta*acceleratorx->neta[j];
        }
        /* jacobian */
        jacx = g4->w[ix] * sqrt(dot(norm, norm));

        /* for each gaussian integration point */
        for (iy = 0; iy < g4->num; iy++)
        {
            xi = g4->xi[iy];
            eta = g4->xi[iy+g4->num];
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                for (s = 0; s < 4; s++)
                    ry[j] += g4->N[iy+s*g4->num]*nodesy[s+4*j];
                norm[j] = acceleratory->n0[j] + xi*acceleratory->nxi[j] + eta*acceleratory->neta[j];
            }
            /* surface normal and jacobian */
            jacy = g4->w[iy] * sqrt(dot(norm, norm));

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr*jacx*jacy;
            *bi += gi*jacx*jacy;
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Compute regular double B integral over a constant TRIA-QUAD element pair */
void dint_triquad_const(const gauss_t *g3,
                        const double *nodesx,
                        const accelerator_t *acceleratorx,
                        const gauss_t *g4,
                        const double *nodesy,
                        const accelerator_t *acceleratory,
                        double k,
                        double *br,
                        double *bi)
{
    int ix, iy, j, s;
    double rx[3], norm[3], jacx;
    double ry[3], jacy, xi, eta;
    double gr, gi;

    *br = *bi = 0.0;

    jacx = sqrt(dot(acceleratorx->n0, acceleratorx->n0));

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each gaussian integration point */
        for (iy = 0; iy < g4->num; iy++)
        {
            xi = g4->xi[iy];
            eta = g4->xi[iy+g4->num];
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 4; s++)
                    ry[j] += g4->N[iy+s*g4->num]*nodesy[s+4*j];
                norm[j] = acceleratory->n0[j] + xi*acceleratory->nxi[j] + eta*acceleratory->neta[j];
            }
            /* surface normal and jacobian */
            jacy = g3->w[ix] * g4->w[iy] * sqrt(dot(norm, norm));

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr*jacy;
            *bi += gi*jacy;
        }
    }

    *br *= jacx;
    *bi *= jacx;
}

/* ------------------------------------------------------------------------ */
/* Compute singular double B integral over a constant TRIA element          */
void dint_tri_const_sing(const gauss_t *g3,
                         const double *nodes,
                         const accelerator_t *accelerator,
                         const gauss_t *g4,
                         double k,
                         double *br,
                         double *bi)
{
    int ix, iy, j, s;
    double rx[3], ry[3], jac;
    double N[3];
    double gr, gi;
    double *wprime, *xiprime, *etaprime;
    int gnum;

    *br = *bi = 0.0;

    wprime = (double *)malloc(sizeof(double)*3*g4->num);
    xiprime = (double *)malloc(sizeof(double)*3*g4->num);
    etaprime = (double *)malloc(sizeof(double)*3*g4->num);

    /* Jacobian over the triangle */
    jac = sqrt(dot(accelerator->n0, accelerator->n0));

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodes[s+3*j];
        }

        sing_quadr_face_tria(g4, g3->xi[ix], g3->xi[ix+g3->num], &gnum, xiprime, etaprime, wprime);

        /* for each gaussian integration point */
        for (iy = 0; iy < gnum; iy++)
        {
            shapefun_tria(xiprime[iy], etaprime[iy], N);

            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 3; s++)
                    ry[j] += N[s]*nodes[s+3*j];
            }

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr * wprime[iy] * g3->w[ix];
            *bi += gi * wprime[iy] * g3->w[ix];
        }
    }

    *br *= jac * jac;
    *bi *= jac * jac;

    free(wprime);
    free(xiprime);
    free(etaprime);
}

/* ------------------------------------------------------------------------ */
/* Compute singular double B integral over a constant QUAD element */
void dint_quad_const_sing(const gauss_t *g4,
                          const double *nodes,
                          const accelerator_t *accelerator,
                          double k,
                          double *br,
                          double *bi)
{
    int ix, iy, j, s, gnum;
    double rx[3], rxi[3], reta[3], norm[3], jacx, jacy;
    double ry[3], xi, eta;
    double gr, gi;
    double N[4];
    double *xiprime, *etaprime, *wprime;

    *br = *bi = 0.0;

    wprime = (double *)malloc(sizeof(double)*4*g4->num);
    xiprime = (double *)malloc(sizeof(double)*4*g4->num);
    etaprime = (double *)malloc(sizeof(double)*4*g4->num);

    /* for each x gaussian integration point */
    for (ix = 0; ix < g4->num; ix++)
    {
        xi = g4->xi[ix];
        eta = g4->xi[ix+g4->num];
        /* computing integration x location and derivatives */
        for (j = 0; j < 3; j++)
        {
            rx[j] = rxi[j] = reta[j] = 0.0;
            for (s = 0; s < 4; s++)
                rx[j] += g4->N[ix+s*g4->num]*nodes[s+4*j];
            norm[j] = accelerator->n0[j] + xi*accelerator->nxi[j] + eta*accelerator->neta[j];
        }
        jacx = g4->w[ix] * sqrt(dot(norm, norm));

        sing_quadr_face_quad(g4, xi, eta, &gnum, xiprime, etaprime, wprime);

        /* for each gaussian integration point */
        for (iy = 0; iy < gnum; iy++)
        {
            shapefun_quad(xiprime[iy], etaprime[iy], N);

            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 4; s++)
                    ry[j] += N[s]*nodes[s+4*j];
                norm[j] = accelerator->n0[j] + xiprime[iy]*accelerator->nxi[j] + etaprime[iy]*accelerator->neta[j];
            }
            /* surface normal and jacobian */
            jacy = wprime[iy] * sqrt(dot(norm, norm));

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr*jacy*jacx;
            *bi += gi*jacy*jacx;
        }
    }

    free(wprime);
    free(xiprime);
    free(etaprime);
}

/* ------------------------------------------------------------------------ */
/* Compute regular double B integral over two linear TRIA elements */
void dint_tritri_lin(const gauss_t *g3,
                     const double *nodesx,
                     const accelerator_t *acceleratorx,
                     const double *nodesy,
                     const accelerator_t *acceleratory,
                     double k,
                     double *br,
                     double *bi)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], jacx, jacy, btempr, btempi;
    double r[3], gr, gi;

    for (s = 0; s < 3*3; s++)
        br[s] = bi[s] = 0.0;

    jacx = sqrt(dot(acceleratorx->n0, acceleratorx->n0));
    jacy = sqrt(dot(acceleratory->n0, acceleratory->n0));

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each y gaussian integration point */
        for (iy = 0; iy < g3->num; iy++)
        {
            /* computing integration y location */
            for (j = 0; j < 3; j++)
            {
                r[j] = -rx[j];
                for (s = 0; s < 3; s++)
                    r[j] += g3->N[iy+s*g3->num]*nodesy[s+3*j];
            }

            green(r, k, &gr, &gi, NULL, NULL, NULL);

            btempr = g3->w[ix] * g3->w[iy] * gr;
            btempi = g3->w[ix] * g3->w[iy] * gi;

            for (qx = 0; qx < 3; qx++)
            {
                for (qy = 0; qy < 3; qy++)
                {
                    br[qx+3*qy] += g3->N[ix+qx*g3->num] * g3->N[iy+qy*g3->num] * btempr;
                    bi[qx+3*qy] += g3->N[ix+qx*g3->num] * g3->N[iy+qy*g3->num] * btempi;
                }
            }
        }
    }

    for (s = 0; s < 3*3; s++)
    {
        br[s] *= jacx*jacy;
        bi[s] *= jacx*jacy;
    }
}

/* ------------------------------------------------------------------------ */
/* Compute regular double B integral over two linear QUAD elements */
void dint_quadquad_lin(const gauss_t *g4,
                       const double *nodesx,
                       const accelerator_t *acceleratorx,
                       const double *nodesy,
                       const accelerator_t *acceleratory,
                       double k,
                       double *br,
                       double *bi)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], ry[3], norm[3], jacx, jacy, xi, eta;
    double gr, gi, btempr, btempi;

    for (s = 0; s < 4*4; s++)
        br[s] = bi[s] = 0.0;

    /* for each x gaussian integration point */
    for (ix = 0; ix < g4->num; ix++)
    {
        xi = g4->xi[ix];
        eta = g4->xi[ix+g4->num];
        /* for each coordinate direction */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 4; s++)
                rx[j] += g4->N[ix+s*g4->num]*nodesx[s+4*j];
            norm[j] = acceleratorx->n0[j] + xi*acceleratorx->nxi[j] + eta*acceleratorx->neta[j];
        }
        /* jacobian */
        jacx = g4->w[ix] * sqrt(dot(norm, norm));

        /* for each gaussian integration point */
        for (iy = 0; iy < g4->num; iy++)
        {
            xi = g4->xi[iy];
            eta = g4->xi[iy+g4->num];
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                for (s = 0; s < 4; s++)
                    ry[j] += g4->N[iy+s*g4->num]*nodesy[s+4*j];
                norm[j] = acceleratory->n0[j] + xi*acceleratory->nxi[j] + eta*acceleratory->neta[j];
            }
            /* surface normal and jacobian */
            jacy = g4->w[iy] * sqrt(dot(norm, norm));

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            btempr = gr*jacx*jacy;
            btempi = gi*jacx*jacy;

            for (qx = 0; qx < 4; qx++)
            {
                for (qy = 0; qy < 4; qy++)
                {
                    br[qx+4*qy] += g4->N[ix+qx*g4->num] * g4->N[iy+qy*g4->num] * btempr;
                    bi[qx+4*qy] += g4->N[ix+qx*g4->num] * g4->N[iy+qy*g4->num] * btempi;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Compute regular double B integral over a linear TRIA-QUAD element pair */
void dint_triquad_lin(const gauss_t *g3,
                      const double *nodesx,
                      const accelerator_t *acceleratorx,
                      const gauss_t *g4,
                      const double *nodesy,
                      const accelerator_t *acceleratory,
                      double k,
                      double *br,
                      double *bi)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], norm[3], jacx;
    double ry[3], jacy, xi, eta;
    double gr, gi;

    for (s = 0; s < 3*4; s++)
        br[s] = bi[s] = 0.0;

    jacx = sqrt(dot(acceleratorx->n0, acceleratorx->n0));

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each gaussian integration point */
        for (iy = 0; iy < g4->num; iy++)
        {
            xi = g4->xi[iy];
            eta = g4->xi[iy+g4->num];
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 4; s++)
                    ry[j] += g4->N[iy+s*g4->num]*nodesy[s+4*j];
                norm[j] = acceleratory->n0[j] + xi*acceleratory->nxi[j] + eta*acceleratory->neta[j];
            }
            /* surface normal and jacobian */
            jacy = g3->w[ix] * g4->w[iy] * sqrt(dot(norm, norm));

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            for (qx = 0; qx < 3; qx++)
            {
                for (qy = 0; qy < 4; qy++)
                {
                    br[qx+3*qy] += g3->N[ix+qx*g3->num] * g4->N[iy+qy*g4->num] * gr*jacy;
                    bi[qx+3*qy] += g3->N[ix+qx*g3->num] * g4->N[iy+qy*g4->num] * gi*jacy;
                }
            }
        }
    }

    for (s = 0; s < 3*4; s++)
    {
        br[s] *= jacx;
        bi[s] *= jacx;
    }
}

/* ------------------------------------------------------------------------ */
/* Compute singular double B integral over a linear TRIA element            */
void dint_tri_lin_sing(const gauss_t *g3,
                       const double *nodes,
                       const accelerator_t *accelerator,
                       const gauss_t *g4,
                       double k,
                       double *br,
                       double *bi)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], ry[3], jac;
    double N[3];
    double gr, gi;
    double *wprime, *xiprime, *etaprime;
    int gnum;

    for (s = 0; s < 3*3; s++)
        br[s] = bi[s] = 0.0;

    wprime = (double *)malloc(sizeof(double)*3*g4->num);
    xiprime = (double *)malloc(sizeof(double)*3*g4->num);
    etaprime = (double *)malloc(sizeof(double)*3*g4->num);

    /* Jacobian over the triangle */
    jac = sqrt(dot(accelerator->n0, accelerator->n0));

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodes[s+3*j];
        }

        sing_quadr_face_tria(g4, g3->xi[ix], g3->xi[ix+g3->num], &gnum, xiprime, etaprime, wprime);

        /* for each gaussian integration point */
        for (iy = 0; iy < gnum; iy++)
        {
            shapefun_tria(xiprime[iy], etaprime[iy], N);

            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 3; s++)
                    ry[j] += N[s]*nodes[s+3*j];
            }

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            for (qx = 0; qx < 3; qx++)
            {
                for (qy = 0; qy < 3; qy++)
                {
                    br[qx+3*qy] += g3->N[ix+qx*g3->num] * g3->w[ix] * N[qy] * wprime[iy] * gr;
                    bi[qx+3*qy] += g3->N[ix+qx*g3->num] * g3->w[ix] * N[qy] * wprime[iy] * gi;
                }
            }
        }
    }

    for (s = 0; s < 3*3; s++)
    {
        br[s] *= jac*jac;
        bi[s] *= jac*jac;
    }

    free(wprime);
    free(xiprime);
    free(etaprime);
}

/* ------------------------------------------------------------------------ */
/* Compute singular double B integral over a linear TRIA element            */
void dint_quad_lin_sing(const gauss_t *g4,
                        const double *nodes,
                        const accelerator_t *accelerator,
                        double k,
                        double *br,
                        double *bi)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], ry[3], jacx, jacy, xi, eta;
    double N[4], norm[3];
    double gr, gi;
    double *wprime, *xiprime, *etaprime;
    int gnum;

    for (s = 0; s < 4*4; s++)
        br[s] = bi[s] = 0.0;

    wprime = (double *)malloc(sizeof(double)*4*g4->num);
    xiprime = (double *)malloc(sizeof(double)*4*g4->num);
    etaprime = (double *)malloc(sizeof(double)*4*g4->num);
	
    /* for each x gaussian integration point */
    for (ix = 0; ix < g4->num; ix++)
    {
        xi = g4->xi[ix];
        eta = g4->xi[ix+g4->num];

        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 4; s++)
                rx[j] += g4->N[ix+s*g4->num]*nodes[s+4*j];
            norm[j] = accelerator->n0[j] + xi * accelerator->nxi[j] + eta * accelerator->neta[j];
        }
        jacx = sqrt(dot(norm, norm)) * g4->w[ix];

        sing_quadr_face_quad(g4, xi, eta, &gnum, xiprime, etaprime, wprime);

        /* for each gaussian integration point */
        for (iy = 0; iy < gnum; iy++)
        {
            xi = xiprime[iy];
            eta = etaprime[iy];

            shapefun_quad(xi, eta, N);
			
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 4; s++)
                    ry[j] += N[s]*nodes[s+4*j];
				
                norm[j] = accelerator->n0[j] + xi * accelerator->nxi[j] + eta * accelerator->neta[j];
            }
            jacy = sqrt(dot(norm, norm)) * wprime[iy];

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            for (qx = 0; qx < 4; qx++)
            {
                for (qy = 0; qy < 4; qy++)
                {
                    br[qx+4*qy] += (g4->N[ix+qx*g4->num] * N[qy]) * gr * (jacx * jacy);
                    bi[qx+4*qy] += (g4->N[ix+qx*g4->num] * N[qy]) * gi * (jacx * jacy);
                }
            }
        }
    }
	
    free(wprime);
    free(xiprime);
    free(etaprime);
}




/* ------------------------------------------------------------------------ */
/* Regular double D integral over two constant TRIA elements                  */
void dintD_tritri_const(const gauss_t *g3,
                       const double *nodesx,
                       const accelerator_t *acceleratorx,
                       const double *nodesy,
                       const accelerator_t *acceleratory,
                       double k,
                       double *br,
                       double *bi)
{
    int ix, iy, j, s;
    double rx[3], r[3], gr, gi;

    *br = *bi = 0.0;

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each y gaussian integration point */
        for (iy = 0; iy < g3->num; iy++)
        {
            /* computing integration y location */
            for (j = 0; j < 3; j++)
            {
                r[j] = -rx[j];
                for (s = 0; s < 3; s++)
                    r[j] += g3->N[iy+s*g3->num]*nodesy[s+3*j];
            }

            green(r, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr*g3->w[ix]*g3->w[iy];
            *bi += gi*g3->w[ix]*g3->w[iy];
        }
    }

	*br *= dot(acceleratorx->n0, acceleratory->n0) * k*k;
	*bi *= dot(acceleratorx->n0, acceleratory->n0) * k*k;
}

/* ------------------------------------------------------------------------ */
/* Compute regular double D integral over two constant QUAD elements        */
void dintD_quadquad_const(const gauss_t *g4,
                         const double *nodesx,
                         const accelerator_t *acceleratorx,
                         const double *nodesy,
                         const accelerator_t *acceleratory,
                         double k,
                         double *br,
                         double *bi)
{
    int ix, iy, j, s;
    double rx[3], ry[3], normx[3], normy[3], xi, eta;
    double gr, gi;

    *br = *bi = 0.0;

    /* for each x gaussian integration point */
    for (ix = 0; ix < g4->num; ix++)
    {
        xi = g4->xi[ix];
        eta = g4->xi[ix+g4->num];
        /* for each coordinate direction */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 4; s++)
                rx[j] += g4->N[ix+s*g4->num]*nodesx[s+4*j];
            normx[j] = acceleratorx->n0[j] + xi*acceleratorx->nxi[j] + eta*acceleratorx->neta[j];
        }

        /* for each gaussian integration point */
        for (iy = 0; iy < g4->num; iy++)
        {
            xi = g4->xi[iy];
            eta = g4->xi[iy+g4->num];
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                for (s = 0; s < 4; s++)
                    ry[j] += g4->N[iy+s*g4->num]*nodesy[s+4*j];
                normy[j] = acceleratory->n0[j] + xi*acceleratory->nxi[j] + eta*acceleratory->neta[j];
            }

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr*g4->w[ix]*g4->w[iy] * dot(normx, normy) * k*k;
            *bi += gi*g4->w[ix]*g4->w[iy] * dot(normx, normy) * k*k;
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Compute regular double D integral over a constant TRIA-QUAD element pair */
void dintD_triquad_const(const gauss_t *g3,
                        const double *nodesx,
                        const accelerator_t *acceleratorx,
                        const gauss_t *g4,
                        const double *nodesy,
                        const accelerator_t *acceleratory,
                        double k,
                        double *br,
                        double *bi)
{
    int ix, iy, j, s;
    double rx[3], normy[3];
    double ry[3], xi, eta;
    double gr, gi;

    *br = *bi = 0.0;

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each gaussian integration point */
        for (iy = 0; iy < g4->num; iy++)
        {
            xi = g4->xi[iy];
            eta = g4->xi[iy+g4->num];
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 4; s++)
                    ry[j] += g4->N[iy+s*g4->num]*nodesy[s+4*j];
                normy[j] = acceleratory->n0[j] + xi*acceleratory->nxi[j] + eta*acceleratory->neta[j];
            }

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr*g3->w[ix] * g4->w[iy] * dot(acceleratorx->n0, normy);
            *bi += gi*g3->w[ix] * g4->w[iy] * dot(acceleratorx->n0, normy);
        }
    }

    *br *= k*k;
    *bi *= k*k;
}

/* ------------------------------------------------------------------------ */
/* Compute singular double D integral over a constant TRIA element          */
void dintD_tri_const_sing(const gauss_t *g3,
                         const double *nodes,
                         const accelerator_t *accelerator,
                         const gauss_t *g4,
                         double k,
                         double *br,
                         double *bi)
{
    int ix, iy, j, s;
    double rx[3], ry[3];
    double N[3];
    double gr, gi;
    double *wprime, *xiprime, *etaprime;
    int gnum;

    *br = *bi = 0.0;

    wprime = (double *)malloc(sizeof(double)*3*g4->num);
    xiprime = (double *)malloc(sizeof(double)*3*g4->num);
    etaprime = (double *)malloc(sizeof(double)*3*g4->num);

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodes[s+3*j];
        }

        sing_quadr_face_tria(g4, g3->xi[ix], g3->xi[ix+g3->num], &gnum, xiprime, etaprime, wprime);

        /* for each gaussian integration point */
        for (iy = 0; iy < gnum; iy++)
        {
            shapefun_tria(xiprime[iy], etaprime[iy], N);

            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 3; s++)
                    ry[j] += N[s]*nodes[s+3*j];
            }

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr * wprime[iy] * g3->w[ix];
            *bi += gi * wprime[iy] * g3->w[ix];
        }
    }

    *br *= dot(accelerator->n0, accelerator->n0) * k*k;
    *bi *= dot(accelerator->n0, accelerator->n0) * k*k;

    free(wprime);
    free(xiprime);
    free(etaprime);
}

/* ------------------------------------------------------------------------ */
/* Compute singular double D integral over a constant QUAD element */
void dintD_quad_const_sing(const gauss_t *g4,
                          const double *nodes,
                          const accelerator_t *accelerator,
                          double k,
                          double *br,
                          double *bi)
{
    int ix, iy, j, s, gnum;
    double rx[3], rxi[3], reta[3], normx[3], normy[3];
    double ry[3], xi, eta;
    double gr, gi;
    double N[4];
    double *xiprime, *etaprime, *wprime;

    *br = *bi = 0.0;

    wprime = (double *)malloc(sizeof(double)*4*g4->num);
    xiprime = (double *)malloc(sizeof(double)*4*g4->num);
    etaprime = (double *)malloc(sizeof(double)*4*g4->num);

    /* for each x gaussian integration point */
    for (ix = 0; ix < g4->num; ix++)
    {
        xi = g4->xi[ix];
        eta = g4->xi[ix+g4->num];
        /* computing integration x location and derivatives */
        for (j = 0; j < 3; j++)
        {
            rx[j] = rxi[j] = reta[j] = 0.0;
            for (s = 0; s < 4; s++)
                rx[j] += g4->N[ix+s*g4->num]*nodes[s+4*j];
            normx[j] = accelerator->n0[j] + xi*accelerator->nxi[j] + eta*accelerator->neta[j];
        }

        sing_quadr_face_quad(g4, xi, eta, &gnum, xiprime, etaprime, wprime);

        /* for each gaussian integration point */
        for (iy = 0; iy < gnum; iy++)
        {
            shapefun_quad(xiprime[iy], etaprime[iy], N);

            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 4; s++)
                    ry[j] += N[s]*nodes[s+4*j];
                normy[j] = accelerator->n0[j] + xiprime[iy]*accelerator->nxi[j] + etaprime[iy]*accelerator->neta[j];
            }

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            *br += gr*wprime[iy]*g4->w[ix] * dot(normx, normy) * k*k;
            *bi += gi*wprime[iy]*g4->w[ix] * dot(normx, normy) * k*k;
        }
    }

    free(wprime);
    free(xiprime);
    free(etaprime);
}






/* ------------------------------------------------------------------------ */
/* Compute regular double D integral over two linear TRIA elements          */
void dintD_tritri_lin(const gauss_t *g3,
                     const double *nodesx,
                     const accelerator_t *acceleratorx,
                     const double *nodesy,
                     const accelerator_t *acceleratory,
                     double k,
                     double *dr,
                     double *di)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], dtempr, dtempi, tempx[3], tempy[3];
    double ry[3], gr, gi;
	double k2nxny, q;

    for (s = 0; s < 3*3; s++)
        dr[s] = di[s] = 0.0;
		
	k2nxny = k*k * dot(acceleratorx->n0, acceleratory->n0);

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each y gaussian integration point */
        for (iy = 0; iy < g3->num; iy++)
        {
            /* computing integration y location */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                for (s = 0; s < 3; s++)
                    ry[j] += g3->N[iy+s*g3->num]*nodesy[s+3*j];
            }

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            dtempr = g3->w[ix] * g3->w[iy] * gr;
            dtempi = g3->w[ix] * g3->w[iy] * gi;

            for (qx = 0; qx < 3; qx++)
            {
				cross(acceleratorx->n0, acceleratorx->gradN + 3*qx, tempx);
                for (qy = 0; qy < 3; qy++)
                {
					cross(acceleratory->n0, acceleratory->gradN + 3*qy, tempy);
					q = k2nxny * g3->N[ix+qx*g3->num] * g3->N[iy+qy*g3->num] - dot(tempx, tempy);
                    dr[qx+3*qy] += q * dtempr;
                    di[qx+3*qy] += q * dtempi;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Compute regular double D integral over two linear TRIA elements          */
void dintD_tritri_lin_fast(const gauss_t *g3,
                     const double *nodesx,
                     const accelerator_t *acceleratorx,
                     const double *nodesy,
                     const accelerator_t *acceleratory,
                     double k,
                     double *dr,
                     double *di)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], dtempr, dtempi, tempx[3], tempy[3];
    double ry[3], gr, gi;
	double q;
	double normx[3], normy[3];
	double jacx, jacy;
	

    for (s = 0; s < 3*3; s++)
        dr[s] = di[s] = 0.0;
		
	jacx = sqrt(dot(acceleratorx->n0, acceleratorx->n0));
	jacy = sqrt(dot(acceleratory->n0, acceleratory->n0));
	for (j = 0; j < 3; j++)
	{
		normx[j] = acceleratorx->n0[j] / jacx;
		normy[j] = acceleratory->n0[j] / jacy;
	}

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each y gaussian integration point */
        for (iy = 0; iy < g3->num; iy++)
        {
            /* computing integration y location */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                for (s = 0; s < 3; s++)
                    ry[j] += g3->N[iy+s*g3->num]*nodesy[s+3*j];
            }

            ddgreen(ry, k, normx, normy, &gr, &gi);
			
			/* mexPrintf("%g", gr); */

            dtempr = g3->w[ix] * g3->w[iy] * gr;
            dtempi = g3->w[ix] * g3->w[iy] * gi;

            for (qx = 0; qx < 3; qx++)
            {
                for (qy = 0; qy < 3; qy++)
                {
					q = g3->N[ix+qx*g3->num] * g3->N[iy+qy*g3->num];
                    dr[qx+3*qy] += q * dtempr;
                    di[qx+3*qy] += q * dtempi;
                }
            }
        }
    }
	
	for (s = 0; s < 3*3; s++)
	{
		dr[s] *= jacx*jacy;
		di[s] *= jacx*jacy;
	}
}



/* ------------------------------------------------------------------------ */
/* Compute regular double B integral over two linear QUAD elements */
void dintD_quadquad_lin(const gauss_t *g4,
                       const double *nodesx,
                       const accelerator_t *acceleratorx,
                       const double *nodesy,
                       const accelerator_t *acceleratory,
                       double k,
                       double *br,
                       double *bi)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], ry[3], norm[3], jacx, jacy, xi, eta;
    double gr, gi, btempr, btempi;

    for (s = 0; s < 4*4; s++)
        br[s] = bi[s] = 0.0;

    /* for each x gaussian integration point */
    for (ix = 0; ix < g4->num; ix++)
    {
        xi = g4->xi[ix];
        eta = g4->xi[ix+g4->num];
        /* for each coordinate direction */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 4; s++)
                rx[j] += g4->N[ix+s*g4->num]*nodesx[s+4*j];
            norm[j] = acceleratorx->n0[j] + xi*acceleratorx->nxi[j] + eta*acceleratorx->neta[j];
        }
        /* jacobian */
        jacx = g4->w[ix] * sqrt(dot(norm, norm));

        /* for each gaussian integration point */
        for (iy = 0; iy < g4->num; iy++)
        {
            xi = g4->xi[iy];
            eta = g4->xi[iy+g4->num];
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                for (s = 0; s < 4; s++)
                    ry[j] += g4->N[iy+s*g4->num]*nodesy[s+4*j];
                norm[j] = acceleratory->n0[j] + xi*acceleratory->nxi[j] + eta*acceleratory->neta[j];
            }
            /* surface normal and jacobian */
            jacy = g4->w[iy] * sqrt(dot(norm, norm));

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            btempr = gr*jacx*jacy;
            btempi = gi*jacx*jacy;

            for (qx = 0; qx < 4; qx++)
            {
                for (qy = 0; qy < 4; qy++)
                {
                    br[qx+4*qy] += g4->N[ix+qx*g4->num] * g4->N[iy+qy*g4->num] * btempr;
                    bi[qx+4*qy] += g4->N[ix+qx*g4->num] * g4->N[iy+qy*g4->num] * btempi;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Compute regular double B integral over a linear TRIA-QUAD element pair */
void dintD_triquad_lin(const gauss_t *g3,
                      const double *nodesx,
                      const accelerator_t *acceleratorx,
                      const gauss_t *g4,
                      const double *nodesy,
                      const accelerator_t *acceleratory,
                      double k,
                      double *br,
                      double *bi)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], norm[3], jacx;
    double ry[3], jacy, xi, eta;
    double gr, gi;

    for (s = 0; s < 3*4; s++)
        br[s] = bi[s] = 0.0;

    jacx = sqrt(dot(acceleratorx->n0, acceleratorx->n0));

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodesx[s+3*j];
        }
        /* for each gaussian integration point */
        for (iy = 0; iy < g4->num; iy++)
        {
            xi = g4->xi[iy];
            eta = g4->xi[iy+g4->num];
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 4; s++)
                    ry[j] += g4->N[iy+s*g4->num]*nodesy[s+4*j];
                norm[j] = acceleratory->n0[j] + xi*acceleratory->nxi[j] + eta*acceleratory->neta[j];
            }
            /* surface normal and jacobian */
            jacy = g3->w[ix] * g4->w[iy] * sqrt(dot(norm, norm));

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            for (qx = 0; qx < 3; qx++)
            {
                for (qy = 0; qy < 4; qy++)
                {
                    br[qx+3*qy] += g3->N[ix+qx*g3->num] * g4->N[iy+qy*g4->num] * gr*jacy;
                    bi[qx+3*qy] += g3->N[ix+qx*g3->num] * g4->N[iy+qy*g4->num] * gi*jacy;
                }
            }
        }
    }

    for (s = 0; s < 3*4; s++)
    {
        br[s] *= jacx;
        bi[s] *= jacx;
    }
}

/* ------------------------------------------------------------------------ */
/* Compute singular double D integral over a linear TRIA element            */
void dintD_tri_lin_sing(const gauss_t *g3,
                       const double *nodes,
                       const accelerator_t *accelerator,
                       const gauss_t *g4,
                       double k,
                       double *dr,
                       double *di)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], ry[3], jac;
    double N[3], tempx[3], tempy[3];
    double gr, gi, dtempr, dtempi;
    double *wprime, *xiprime, *etaprime, q, k2nxny;
    int gnum;

    for (s = 0; s < 3*3; s++)
        dr[s] = di[s] = 0.0;

    wprime = (double *)malloc(sizeof(double)*3*g4->num);
    xiprime = (double *)malloc(sizeof(double)*3*g4->num);
    etaprime = (double *)malloc(sizeof(double)*3*g4->num);
	
	k2nxny = k*k * dot(accelerator->n0, accelerator->n0);

    /* for each x gaussian integration point */
    for (ix = 0; ix < g3->num; ix++)
    {
        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 3; s++)
                rx[j] += g3->N[ix+s*g3->num]*nodes[s+3*j];
        }

        sing_quadr_face_tria(g4, g3->xi[ix], g3->xi[ix+g3->num], &gnum, xiprime, etaprime, wprime);

        /* for each gaussian integration point */
        for (iy = 0; iy < gnum; iy++)
        {
            shapefun_tria(xiprime[iy], etaprime[iy], N);

            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 3; s++)
                    ry[j] += N[s]*nodes[s+3*j];
            }

            green(ry, k, &gr, &gi, NULL, NULL, NULL);
			
            dtempr = g3->w[ix] * wprime[iy] * gr;
            dtempi = g3->w[ix] * wprime[iy] * gi;

            for (qx = 0; qx < 3; qx++)
            {
				cross(accelerator->n0, accelerator->gradN + 3*qx, tempx);
                for (qy = 0; qy < 3; qy++)
                {
					cross(accelerator->n0, accelerator->gradN + 3*qy, tempy);
					q = (k2nxny * g3->N[ix+qx*g3->num] * N[qy] - dot(tempx, tempy));
                    dr[qx+3*qy] += q * dtempr;
                    di[qx+3*qy] += q * dtempi;
                }
            }
        }
    }

    free(wprime);
    free(xiprime);
    free(etaprime);
}

/* ------------------------------------------------------------------------ */
/* Compute singular double B integral over a linear TRIA element            */
void dintD_quad_lin_sing(const gauss_t *g4,
                        const double *nodes,
                        const accelerator_t *accelerator,
                        double k,
                        double *br,
                        double *bi)
{
    int ix, iy, j, s, qx, qy;
    double rx[3], ry[3], jacx, jacy, xi, eta;
    double N[4], norm[3];
    double gr, gi;
    double *wprime, *xiprime, *etaprime;
    int gnum;

    for (s = 0; s < 4*4; s++)
        br[s] = bi[s] = 0.0;

    wprime = (double *)malloc(sizeof(double)*4*g4->num);
    xiprime = (double *)malloc(sizeof(double)*4*g4->num);
    etaprime = (double *)malloc(sizeof(double)*4*g4->num);
	
    /* for each x gaussian integration point */
    for (ix = 0; ix < g4->num; ix++)
    {
        xi = g4->xi[ix];
        eta = g4->xi[ix+g4->num];

        /* computing integration x location */
        for (j = 0; j < 3; j++)
        {
            rx[j] = 0.0;
            for (s = 0; s < 4; s++)
                rx[j] += g4->N[ix+s*g4->num]*nodes[s+4*j];
            norm[j] = accelerator->n0[j] + xi * accelerator->nxi[j] + eta * accelerator->neta[j];
        }
        jacx = sqrt(dot(norm, norm)) * g4->w[ix];

        sing_quadr_face_quad(g4, xi, eta, &gnum, xiprime, etaprime, wprime);

        /* for each gaussian integration point */
        for (iy = 0; iy < gnum; iy++)
        {
            xi = xiprime[iy];
            eta = etaprime[iy];

            shapefun_quad(xi, eta, N);
			
            /* for each coordinate direction */
            for (j = 0; j < 3; j++)
            {
                ry[j] = -rx[j];
                /* computing integration location and its derivatives */
                for (s = 0; s < 4; s++)
                    ry[j] += N[s]*nodes[s+4*j];
				
                norm[j] = accelerator->n0[j] + xi * accelerator->nxi[j] + eta * accelerator->neta[j];
            }
            jacy = sqrt(dot(norm, norm)) * wprime[iy];

            green(ry, k, &gr, &gi, NULL, NULL, NULL);

            for (qx = 0; qx < 4; qx++)
            {
                for (qy = 0; qy < 4; qy++)
                {
                    br[qx+4*qy] += (g4->N[ix+qx*g4->num] * N[qy]) * gr * (jacx * jacy);
                    bi[qx+4*qy] += (g4->N[ix+qx*g4->num] * N[qy]) * gi * (jacx * jacy);
                }
            }
        }
    }
	
    free(wprime);
    free(xiprime);
    free(etaprime);
}
