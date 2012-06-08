#ifndef FMBEM_H
#define FMBEM_H

#include <stdint.h>


void leafGq(int nnod,
            const double *r,
            const double *qr,
            const double *qi,
            const int32_t *father,
            int nclus,
            int ns,
            const double *s,
            double k,
            double *Fr,
            double *Fi);

void leafHp(int nnod,
            const double *r,
            const double *n,
            const double *pr,
            const double *pi,
            const int32_t *father,
            int nclus,
            int ns,
            const double *s,
            double k,
            double *Fr,
            double *Fi);

void recover(int nnod,
             const double *r,
             const double *Nr,
             const double *Ni,
             const int32_t *father,
             int ns,
             const double *s,
             const double *w,
             double k,
             double *pr,
             double *pi);

void upward(int nnod,
            const double *r,
            const double *qr,
            const double *qi,
            const int32_t *father,
            int nclus,
            int ns,
            const double *s,
            double k,
            double *Fr,
            double *Fi);

void translate(int nclus,
               int ns, const double *Fr,
               const double *Fi,
               int nil,
               const int32_t *I,
               const int32_t *D,
               const int32_t *P,
               const int32_t *Perm,
               int nm,
               const double *Mr,
               const double *Mi,
               double *Nr,
               double *Ni);

void downward(int nnod,
              const double *r,
              const double *qr,
              const double *qi,
              const int32_t *father,
              int ns,
              const double *s,
              double k,
              double *Nr,
              double *Ni);
#endif
