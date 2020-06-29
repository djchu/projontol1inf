#ifndef SemismoothNewton_H
#define SemismoothNewton_H

void print_Data(double *m, int nRows,int nCols);
double norm_l1infinity(double *m, int nRows, int nCols);

void SemismoothNewton(double B[], double C, double A[],
                        int nRows, int nCols);
#endif // PROJL1INFNEWTON_H

