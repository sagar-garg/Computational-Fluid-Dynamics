#include "sor.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  int **flag,
  int number_fluid
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      if(flag[i][j]==1){
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
    }
  }

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
        if(flag[i][j]==1){
        rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
        }
    }
  }
  rloc = rloc/(number_fluid);//Not to divide with (imax*jmax)
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;

/* set boundary values */
//How it is encounting if there is a inflow at one of the boundary>?
  for(i = 1; i <= imax; i++) {
    P[i][0] = P[i][1];//Lower boundary 
    P[i][jmax+1] = P[i][jmax];//Top boundary
  }
  for(j = 1; j <= jmax; j++) {
    P[0][j] = P[1][j];//Left boundary
    P[imax+1][j] = P[imax][j];//Right boundary
  }


for(int i = 1; i <=imax; i++){
    for (int j = 1; j <=jmax; j++)
    {
        switch(flag[i][j]){
            case 66: case 96:   P[i][j] = P[i][j+1];//BN
                break;
            case 258: case 288: P[i][j] = P[i-1][j];//BW
                break;
            case 130: case 160: P[i][j] = P[i][j-1];//BS
                break;
            case 514: case 544: P[i][j] = P[i+1][j];//BO
                break;
            case 578: case 608: P[i][j] = (P[i][j+1] + P[i+1][j])/2;//BNO
                break;
            case 322: case 352: P[i][j] = (P[i][j+1] + P[i-1][j])/2;//BNW
                break;
            case 386: case 416: P[i][j] = (P[i-1][j] + P[i][j-1])/2;//BSW
                break;
            case 642: case 670: P[i][j] = (P[i][j-1] + P[i+1][j])/2;//BSO
                break;
    }
}
  
}
}
