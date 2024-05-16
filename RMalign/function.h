#ifndef function_H_H
#define function_H_H
double dist(double x[3], double y[3]);
double dot(double *a, double *b);
void transform(double t[3], double u[3][3], double *x, double *x1);
void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3]);
#endif
