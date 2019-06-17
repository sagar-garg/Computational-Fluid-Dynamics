//  boundary_val.c

#include "boundary_val.h"
#include <stdio.h>

void boundaryvalues(
    int imax,
    int jmax,
    double **U,
    double **V)
{
    for(int i = 1; i <= imax; i++){
        V[i][0] = 0;
        V[i][jmax] = 0;
        U[i][0] = -U[i][1];
        U[i][jmax+1] = 2.0 - U[i][jmax];
        
    }
    
    for(int j = 1; j <= jmax; j++){
        U[0][j] = 0;
        U[imax][j] =0;
        V[0][j] = -V[1][j];
        V[imax+1][j] = -V[imax][j];
    }

}
