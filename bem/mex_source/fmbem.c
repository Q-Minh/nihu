/* fmbem.c  implementation of the fast multipole bem C routines */
#include "fmbem.h"

#include "math.h"
#include "vector.h"

/* Compute the leaf level far field signatures of the product Gq */
void leafGq(int nnod,               /* number of nodes */
            const double *r,        /* nodal coordinates */
            const double *qr,       /* real part of nodal velocities */
            const double *qi,       /* imaginary part of nodal velocities */
            const int32_t *father,  /* father cluster indices of nodes */
            int nclus,              /* number of clusters on leaf level */
            int ns,                 /* number of points on the unit sphere */
            const double *s,        /* points on the unit sphere */
            double k,               /* real wave number */
            double *Fr,             /* real part of far field signarures */
            double *Fi)             /* imaginary part of far field signarures */
{
    /* Clear output matrix */
    int j, p;
    for (j = 0; j < nclus; j++)
        for (p = 0; p < ns; p++)
            Fr[j*ns+p] = Fi[j*ns+p] = 0.0;

    /* For each source node */
    for (j = 0; j < nnod; j++)
    {
        /* father cluster index */
        int32_t f = father[j];
        /* for each direction (s) */
        for (p = 0; p < ns; p++)
        {
            /* F = q * exp(-i*k*d*s) */
            /* phi = d*s */
            double phi = k * dot(&r[j*3], &s[p*3]);
            double cphi = cos(phi);
            double sphi = sin(phi);
            Fr[f*ns+p] += qr[j]*cphi+qi[j]*sphi;
            Fi[f*ns+p] += qi[j]*cphi-qr[j]*sphi;
        }
    }
}

/* Compute the leaf level far field signatures of the product Hp */
void leafHp(int nnod,
            const double *r,        /* nodal coordinates */
            const double *n,        /* nodal normals */
            const double *pr,       /* real part of nodal pressures */
            const double *pi,       /* imaginary part of nodal pressures */
            const int32_t *father,  /* father cluster indices of nodes */
            int nclus,              /* number of clusters on leaf level */
            int ns,                 /* number of points on the unit sphere */
            const double *s,        /* points on the unit sphere */
            double k,               /* real wave number */
            double *Fr,             /* real part of far field signarures */
            double *Fi)             /* imaginary part of far field signarures */
{
    /* Clear output matrix */
    int j, p;
    for (j = 0; j < nclus; j++)
        for (p = 0; p < ns; p++)
            Fr[j*ns+p] = Fi[j*ns+p] = 0.0;

    /* For each source node */
    for (j = 0; j < nnod; j++)
    {
        int i;
        double kdvec[3];
        /* father cluster index */
       	int32_t f = father[j];
        /* distance vector (kd) */
        for (i = 0; i < 3; i++)
            kdvec[i] = k*r[j*3+i];
        /* for each	direction (s) */
        for (p = 0; p < ns; p++)
        {
            /* phi = d*s */
            double phi = dot(kdvec, &s[p*3]);
            double cphi = cos(phi);
            double sphi = sin(phi);
            /* F = i*k* q * exp(-i*k*d*s) * (n*s) */
            /* k*n*s */
            double kns = k*dot(&n[j*3], &s[p*3]);
            Fr[f*ns+p] -= kns*(pi[j]*cphi-pr[j]*sphi);
            Fi[f*ns+p] += kns*(pr[j]*cphi+pi[j]*sphi);
        }
    }
}

void recover(int nnod,              /* number of nodes */
             const double *r,       /* nodal distances */
             const double *Nr,      /* real part of near field signatures */
             const double *Ni,      /* imaginary part of near field signatures */
             const int32_t *father, /* father cluster indices */
             int ns,                /* number of points on the unit sphere */
             const double *s,       /* points on the unit sphere */
             const double *w,       /* quadrature weights on the unit sphere */
             double k,              /* wave number */
             double *pr,            /* real part of pressure */
             double *pi)            /* imaginary part of pressure */
{
    /* for each receiver node */
    int j;
    for (j = 0; j < nnod; j++)
    {
        double kdvec[3];
        int i, p;
        /* father cluster of receiver */
        int32_t f = father[j];
        /* numerical integration initialization */
        pr[j] = pi[j] = 0.0;
        /* distance vector ( k*(r-R) ) */
        for (i = 0; i < 3; i++)
            kdvec[i] = -k*(r[j*3+i]);
        /* for each direction */
        for (p = 0; p < ns; p++)
        {
            /* phi = k*d*s */
            double phi = dot(kdvec, &s[p*3]);
            double cphi = cos(phi);
            double sphi = sin(phi);
            /* p += N*exp(-i*phi)*w */
            pr[j] += (Nr[f*ns+p]*cphi+Ni[f*ns+p]*sphi) * w[p];
            pi[j] += (Ni[f*ns+p]*cphi-Nr[f*ns+p]*sphi) * w[p];
        }
    }
}

void upward(int nnod,               /* number of nodes */
            const double *r,        /* nodal distances */
            const double *qr,       /* real part of signatures */
            const double *qi,       /* imaginary part of signatures */
            const int32_t *father,  /* father indices */
            int nclus,              /* number of clusters */
            int ns,                 /* number of points on the unit sphere */
            const double *s,        /* points on the unit sphere */
            double k,               /* wave number */
            double *Fr,             /* real part of far field signatures */
            double *Fi)             /* imaginary part of far field signatures */
{
    int j, p;
    for (j = 0; j < nclus; j++)
        for (p = 0; p < ns; p++)
            Fr[j*ns+p] = Fi[j*ns+p] = 0.0;

    for (j = 0; j < nnod; j++)
    {
        int32_t f = father[j];
        double kdvec[3];
        int i;
        for (i = 0; i < 3; i++)
            kdvec[i] = k*r[j*3+i];
        for (p = 0; p < ns; p++)
        {
            double phi = dot(kdvec, &s[3*p]);
            double cphi = cos(phi);
            double sphi = sin(phi);
            Fr[f*ns+p] += qr[j*ns+p]*cphi+qi[j*ns+p]*sphi;
            Fi[f*ns+p] += qi[j*ns+p]*cphi-qr[j*ns+p]*sphi;
        }
    }
}

void translate(int nclus,           /* number of clusters */
               int ns,              /* number of points on the unit sphere */
               const double *Fr,    /* real part of far field signatures */
               const double *Fi,    /* imaginary part of far field signatures */
               int nil,             /* max. length of interaction list */
               const int32_t *I,    /* interaction list matrix */
               const int32_t *D,    /* distance matrix */
               const int32_t *P,    /* permutation index matrix */
               const int32_t *Perm, /* permutation matrix */
               int nm,              /* number of translation operators */
               const double *Mr,    /* real part of translation operators */
               const double *Mi,    /* imaginary part of translation operators */
               double *Nr,          /* real part of near field signatures */
               double *Ni)          /* imaginary part of near field signatures */
{
    /* Clear output matrix */
    int i, p;
    for (i = 0; i < nclus; i++)
        for (p = 0; p < ns; p++)
            Nr[i*ns+p] = Ni[i*ns+p] = 0.0;

    /* for each source cluster */
    for (i = 0; i < nclus; i++)
    {
        /* for each receiver cluster in the interaction list */
        int j;
        for (j = 0; j < nil; j++)
        {
            /* receiver cluster index */
            int32_t Iij = I[i*nil+j];
            if (Iij >= 0)
            {
                /* distance index */
                int32_t Dij = D[i*nil+j];
                /* permutation index */
                int32_t Pij = P[i*nil+j];
                /* for each direction (unit sphere) */
                for (p = 0; p < ns; p++)
                {
                    /* permutated index */
                    int32_t pp = Perm[Pij*ns+p];
                    /* source */
                    double fr = Fr[i*ns+p];
                    double fi = Fi[i*ns+p];
                    /* translation */
                    double mr = Mr[Dij*ns+pp];
                    double mi = Mi[Dij*ns+pp];
                    /* product */
                    Nr[Iij*ns+p] += fr*mr-fi*mi;
                    Ni[Iij*ns+p] += fr*mi+fi*mr;
                }
            }
        }
    }
}

void downward(int nnod,             /* number of nodes */
              const double *r,      /* nodal distances */
              const double *qr,     /* real part of father signatures */
              const double *qi,     /* imaginary part of father signatures */
              const int32_t *father,/* father cluster indices */
              int ns,               /* number of points on the unit sphere */
              const double *s,      /* points on the unit sphere */
              double k,             /* wave number */
              double *Nr,           /* real part of near field signatures */
              double *Ni)           /* imaginary part of near field signatures */
{
    /* for each receiver cluster */
    int j;
    for (j = 0; j < nnod; j++)
    {
        int i, p;
        /* index of father source cluster */
        int32_t f = father[j];
        /* distance vector (k*(R-r)) */
        double kdvec[3];
        for (i = 0; i < 3; i++)
            kdvec[i] = -k*r[j*3+i];
        /* for each direction (unit sphere) */
        for (p = 0; p < ns; p++)
        {
            /* phi = kds */
            double phi = dot(kdvec, &s[3*p]);
            double cphi = cos(phi);
            double sphi = sin(phi);
            /* N = q*exp(-kds) */
            Nr[j*ns+p] = qr[f*ns+p]*cphi+qi[f*ns+p]*sphi;
            Ni[j*ns+p] = qi[f*ns+p]*cphi-qr[f*ns+p]*sphi;
        }
    }
}
