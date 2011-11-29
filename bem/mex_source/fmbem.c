#include "fmbem.h"
#include "math.h"
#include "vector.h"

void leafGq(int nnod,
            const double *r,
            const double *qr,
            const double *qi,
            const double *father,
            int nclus,
            int ns,
            const double *s,
            double k,
            double *Fr,
            double *Fi)
{
    int j, p, f;
    double phi, cphi, sphi;

    /* Clear output matrix */
    for (j = 0; j < nclus; j++)
        for (p = 0; p < ns; p++)
            Fr[j*ns+p] = Fi[j*ns+p] = 0.0;

    /* For each source node */
    for (j = 0; j < nnod; j++)
    {
        /* father cluster index */
        f = (int)father[j];
        /* for each	direction (s) */
        for (p = 0; p < ns; p++)
        {
            /* F = q * exp(-i*k*d*s) */
            /* phi = d*s */
            phi = k * dot(&r[j*3], &s[p*3]);
            cphi = cos(phi);
            sphi = sin(phi);
            Fr[f*ns+p] += qr[j]*cphi+qi[j]*sphi;
            Fi[f*ns+p] += qi[j]*cphi-qr[j]*sphi;
        }
    }
}

void leafHp(int nnod,
            const double *r,
            const double *n,
            const double *pr,
            const double *pi,
            const double *father,
            int nclus,
            int ns,
            const double *s,
            double k,
            double *Fr,
            double *Fi)
{
    int i, j, p, f;
    double kdvec[3], phi, cphi, sphi, kns;

    /* Clear output matrix */
    for (j = 0; j < nclus; j++)
        for (p = 0; p < ns; p++)
            Fr[j*ns+p] = Fi[j*ns+p] = 0.0;

    /* For each source node */
    for (j = 0; j < nnod; j++)
    {
        /* father cluster index */
        f = (int)father[j];
        /* distance vector (kd) */
        for (i = 0; i < 3; i++)
            kdvec[i] = k*r[j*3+i];
        /* for each	direction (s) */
        for (p = 0; p < ns; p++)
        {
            /* F = i*k* q * exp(-i*k*d*s) * (n*s) */
            /* k*n*s */
            kns = k*dot(&n[j*3], &s[p*3]);
            /* phi = d*s */
            phi = dot(kdvec, &s[p*3]);
            cphi = cos(phi);
            sphi = sin(phi);
            Fr[f*ns+p] -= kns*(pi[j]*cphi-pr[j]*sphi);
            Fi[f*ns+p] += kns*(pr[j]*cphi+pi[j]*sphi);
        }
    }
}

void recover(int nnod,
             const double *r,
             const double *Nr,
             const double *Ni,
             const double *father,
             int ns,
             const double *s,
             const double *w,
             double k,
             double *pr,
             double *pi)
{
    int i, j, p, f;
    double kdvec[3], phi, cphi, sphi;

    /* for each receiver node */
    for (j = 0; j < nnod; j++)
    {
        /* numerical integration initialization */
        pr[j] = pi[j] = 0.0;
        /* father cluster of receiver */
        f = (int)father[j];
        /* distance vector ( k*(r-R) ) */
        for (i = 0; i < 3; i++)
            kdvec[i] = -k*(r[j*3+i]);
        /* for each direction */
        for (p = 0; p < ns; p++)
        {
            /* phi = k*d*s */
            phi = dot(kdvec, &s[p*3]);
            cphi = cos(phi);
            sphi = sin(phi);
            /* p += N*exp(-phi)*w */
            pr[j] += (Nr[f*ns+p]*cphi+Ni[f*ns+p]*sphi) * w[p];
            pi[j] += (Ni[f*ns+p]*cphi-Nr[f*ns+p]*sphi) * w[p];
        }
    }
}

void upward(int nnod,
            const double *r,
            const double *qr,
            const double *qi,
            const double *father,
            int nclus,
            int ns,
            const double *s,
            double k,
            double *Fr,
            double *Fi)
{
    int i, j, p, f;
    double kdvec[3], phi, cphi, sphi;

    for (j = 0; j < nclus; j++)
        for (p = 0; p < ns; p++)
            Fr[j*ns+p] = Fi[j*ns+p] = 0.0;

    for (j = 0; j < nnod; j++)
    {
        f = (int)father[j];
        for (i = 0; i < 3; i++)
            kdvec[i] = k*r[j*3+i];
        for (p = 0; p < ns; p++)
        {
            phi = dot(kdvec, &s[3*p]);
            cphi = cos(phi);
            sphi = sin(phi);
            Fr[f*ns+p] += qr[j*ns+p]*cphi+qi[j*ns+p]*sphi;
            Fi[f*ns+p] += qi[j*ns+p]*cphi-qr[j*ns+p]*sphi;
        }
    }
}

void translate(int nclus,
               int ns,
               const double *Fr,
               const double *Fi,
               int nil,
               const double *I,
               const double *D,
               const double *P,
               const double *Perm,
               int nm,
               const double *Mr,
               const double *Mi,
               double *Nr,
               double *Ni)
{
    int i, j, p, pp, Iij, Dij, Pij;
    double fr, fi, mr, mi;

    /* Clear output matrix */
    for (i = 0; i < nclus; i++)
        for (p = 0; p < ns; p++)
            Nr[i*ns+p] = Ni[i*ns+p] = 0.0;

    /* for each source cluster */
    for (i = 0; i < nclus; i++)
    {
        /* for each receiver cluster in the interaction list */
        for (j = 0; j < nil; j++)
        {
            /* receiver cluster index */
            Iij = (int)I[i*nil+j];
            if (Iij >= 0)
            {
                /* distance index */
                Dij = (int)D[i*nil+j];
                /* permutation index */
                Pij = (int)P[i*nil+j];
                /* for each direction (unit sphere) */
                for (p = 0; p < ns; p++)
                {
                    /* permutated index */
                    pp = (int)Perm[Pij*ns+p];
                    /* source */
                    fr = Fr[i*ns+p];
                    fi = Fi[i*ns+p];
                    /* translation */
                    mr = Mr[Dij*ns+pp];
                    mi = Mi[Dij*ns+pp];
                    /* product */
                    Nr[Iij*ns+p] += fr*mr-fi*mi;
                    Ni[Iij*ns+p] += fr*mi+fi*mr;
                }
            }
        }
    }
}

void downward(int nnod,
              const double *r,
              const double *qr,
              const double *qi,
              const double *father,
              int ns,
              const double *s,
              double k,
              double *Nr,
              double *Ni)
{
    int i, j, p, f;
    double kdvec[3], phi, cphi, sphi;

    /* for each receiver cluster */
    for (j = 0; j < nnod; j++)
    {
        /* index of father source cluster */
        f = (int)father[j];
        /* distance vector (k*(R-r)) */
        for (i = 0; i < 3; i++)
            kdvec[i] = -k*r[j*3+i];
        /* for each direction (unit sphere) */
        for (p = 0; p < ns; p++)
        {
            /* phi = kds */
            phi = dot(kdvec, &s[3*p]);
            cphi = cos(phi);
            sphi = sin(phi);
            /* N = q*exp(-kds) */
            Nr[j*ns+p] = qr[f*ns+p]*cphi+qi[f*ns+p]*sphi;
            Ni[j*ns+p] = qi[f*ns+p]*cphi-qr[f*ns+p]*sphi;
        }
    }
}
