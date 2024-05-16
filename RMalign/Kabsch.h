#ifndef Kasch_H_H
#define Kasch_H_H
#include<math.h>
bool Kabsch(double **x, 
            double **y, 
            int n, 
            int mode,
            double *rms,
            double t[3],
            double u[3][3]
            );
#endif
