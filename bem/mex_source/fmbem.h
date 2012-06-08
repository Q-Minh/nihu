/* fmbem.h  header file for the fast multipole bem C routines */
#ifndef FMBEM_H
#define FMBEM_H

#include <stdint.h> /* for the definition of int32_t */

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
            double *Fi);            /* imaginary part of far field signarures */

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
            double *Fi);            /* imaginary part of far field signarures */

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
             double *pi);           /* imaginary part of pressure */

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
            double *Fi);            /* imaginary part of far field signatures */

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
               double *Ni);         /* imaginary part of near field signatures */

void downward(int nnod,             /* number of nodes */
              const double *r,      /* nodal distances */
              const double *qr,     /* real part of father signatures */
              const double *qi,     /* imaginary part of father signatures */
              const int32_t *father,/* father cluster indices */
              int ns,               /* number of points on the unit sphere */
              const double *s,      /* points on the unit sphere */
              double k,             /* wave number */
              double *Nr,           /* real part of near field signatures */
              double *Ni);          /* imaginary part of near field signatures */
#endif /* ifndef FMBEM_H */
