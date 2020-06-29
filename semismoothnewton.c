#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "semismoothnewton.h"

/*
 * compute the value of F(x)
 */
void fmutheta(double Fval[], double grad_mu[], double B[],
              double mu[], int reduced_inx[],
              int reduced_Rows, int nRows, int nCols,
              double theta,
              double C){
    int irow, index;
    double mu_i, sum_mu = 0.0;
    for(int i=0; i<reduced_Rows; i++){
        irow = reduced_inx[i];
        mu_i = mu[irow];
        sum_mu += mu_i;

        grad_mu[irow] = 0.0;
        Fval[irow] = -theta;
        double maxAmu0 = 0.0;
        for(int j=0; j<nCols; j++){
            index = irow*nCols + j;
            maxAmu0 = B[index]-mu_i;

            if(maxAmu0 >= 0.0){
                Fval[irow] += maxAmu0;
                grad_mu[irow] -= 1.0;
            }
        }
    }
    Fval[nRows] = -sum_mu + C;
}

void print_Data(double *m, int nRows,int nCols) {
    int i,j,index;
    double value;

    for(i=0;i<nRows;i++) {
        for (j=0;j<nCols;j++) {
            index = i*nCols + j;
            value = m[index];
            printf("%f ",value);
        }
        printf("\n");
    }
}

double norm_l1infinity(double *m, int nRows, int nCols) {
    double norm=0;
    int i, j, index;
    for(i=0;i<nRows;i++) {
        double maxRow=0;
        for (j=0;j<nCols;j++) {
            index=i*nCols + j;
            double value=fabs(m[index]);
            if(value>maxRow)	{
                maxRow=value;
            }
            printf(" %.4f ", value);
        }
        fflush(stdout);
        norm= norm + maxRow;
    }
    return norm;
}

/*
 * Semismooth Newton Algorithm for the L_{1,inf}-norm projection
 */
void mySemismoothNewton(double B[], double C, double A[],
                        int nRows, int nCols){

    double theta = 0.0, L1Infnorm = 0.0, d_norm = 0.0;
    int i, j, iter=0, maxIter = 30;
    double *mu = (double *) malloc(sizeof(double)*nRows);
    double *thetaBound = (double *) malloc(sizeof(double)*nRows);
    
    int *reduced_inx = (int *)malloc(sizeof(int)*nRows);
    int reduced_Rows;
    double *Fval = (double *)malloc(sizeof(double)*(nRows+1)); // F(\mu, \theta)
    double *grad_mu = (double *)malloc(sizeof(double)*nRows); // n_i
    double *d = (double *)malloc(sizeof(double)*(nRows+1)); // search direction
    double *z = (double *)malloc(sizeof(double)*(nRows+1)); // LU decomposition

    // Initialize the variables
    for(i=0; i<nRows; i++){
        Fval[i] = 0.0;
        grad_mu[i] = 0.0;
        d[i] = 0.0;
        z[i] = 0.0;
    }
    Fval[nRows] = 0.0;
    d[nRows] = 0.0;
    z[nRows] = 0.0;

//    printf("A: ");
//    for(i=0; i<nRows*nCols; i++)
//        printf("%f ", A[i]);
//    printf("\n");

    // Step 1: set the initial value to mu and
    //         compute the L_{1,inf} norm of matrix A
    for( i=0; i<nRows; i++){
        double maxAi = fabs(A[i*nCols]);
        double sumAi = 0.0;
//         maxAi=fabs(A[i*nCols]);
//         printf("MAX of Row %d: ", i);
        for( j=0; j<nCols; j++){
            int index = i*nCols + j;
            B[index] = A[index];
            double absBinx = fabs(B[index]);

            if(maxAi<absBinx)
                maxAi = absBinx;
            sumAi = sumAi + absBinx; 
        }
//        printf("%f\n", maxAi);
        mu[i] = maxAi;  //initialize mu
        thetaBound[i] = sumAi;
        L1Infnorm += maxAi;
    }

    if(L1Infnorm<=C){
        free(mu);       free(thetaBound);
        free(reduced_inx);
        free(Fval);     free(grad_mu);
        free(d);        free(z);
        return; //Return matrix B
    }

    for( i=0; i<nRows*nCols; i++){
        if(B[i]<0)
            B[i] = fabs(B[i]);
    }

    // Step 2: Newton Algorithm to solve the system
    //Initialize the row index
    for(j=0, i=0; i<nRows; i++){    
        reduced_inx[j] = i;
        j++;
    }
    reduced_Rows = j;

    //printf("Semismooth Newton Algorithm begins ...\n");
    while(iter<maxIter){

        // compute F(mu, theta)
        fmutheta(Fval, grad_mu, B, mu, reduced_inx,
                 reduced_Rows, nRows, nCols, theta, C);
        /*
        printf("Fval: ");
        for(int i=0; i<=nRows; i++)
            printf("%f; ", Fval[i]);
        printf("\n");
         */

        // compute d via LU decomposition instead of inverse matrix of Jacobian
        int irow;
        double sumz = 0.0;
        double sum_gradmu = 0.0;
        for( i=0; i<reduced_Rows; i++){
            irow = reduced_inx[i];
            z[irow] = -Fval[irow];
            sumz += z[irow]/grad_mu[irow];
            sum_gradmu -= 1.0/grad_mu[irow];
        }
        z[nRows] = -Fval[nRows] + sumz;

        d[nRows] = z[nRows]/sum_gradmu;
        d_norm = d[nRows]*d[nRows];
        for( i=0; i<reduced_Rows; i++){
            irow = reduced_inx[i];
            d[irow] = (z[irow] + d[nRows])/grad_mu[irow];
            d_norm += d[irow]*d[irow];
        }
        d_norm = sqrt(d_norm);
        //printf("iter: %d, norm(d): %f\n", iter, d_norm);
        iter++;

        if(d_norm<1.0e-10){
           break;
        }

        for( i=0; i<reduced_Rows; i++){
            irow = reduced_inx[i];
            mu[irow] = mu[irow] + d[irow];
            //if(mu[irow]<0)
            //    mu[irow] = 0;
            //if(mu[irow]>maxAi[irow])
            //    mu[irow] = maxAi[irow];
            
            //printf(" mu[%d]: %f ", irow, mu[irow]);        
        }
        theta = theta + d[nRows];
        //printf(" theta: %f\n", theta);    
        
        for(j=0, i=0; i<nRows; i++){
            if(theta <= thetaBound[i]){
                reduced_inx[j] = i;
                j++;
            }
        }
        reduced_Rows = j;
    }

    // Step 3: Output the projection B
    for( i=0; i<nRows; i++){
        int index;
        //printf("Mu[%d]=%f, ", i, mu[i]);
        if(mu[i]<0)
            mu[i] = 0;
        
        for( j=0; j<nCols; j++){
            index = i*nCols + j;
            if(A[index]<0)
                B[index] = -fmin(B[index], mu[i]);
            else
                B[index] = fmin(B[index], mu[i]);
//            printf("B[%d]=%f, ", index, B[index]);
        }
//        printf("\n");
    }
//    printf("\n");

    free(mu);      free(thetaBound);
    free(reduced_inx);
    free(Fval);     free(grad_mu);
    free(d);        free(z);

}
