//
//  uvp.c
//  CFD_Lab_test1
//
//  Created by Sagar Garg on 06/05/19.
//  Copyright © 2019 Sagar Garg. All rights reserved.
//

#include <stdio.h>
#include <math.h>

#include "helper.h"
// #include "init.h"
#include "uvp.h"



// declaring functions used for discretization
double du_dx( double uim1j, double uij, double dx);
double d2u_dx2( double uim1j, double uij, double uip1j, double dx);

double duv_dx (double uim1j, double uim1jp1, double uij, double uijp1, double vim1j, double vij, double vip1j, double dx, double alpha);
double du2_dx (double uim1j, double uij, double uip1j, double dx, double alpha);

//helper function for calculate_dt. find maximum value in a matrix
double max_val_U(int imax,int jmax,double** U);
double max_val_V(int imax,int jmax,double** V);


void calculate_fg(
                  double Re,
                  double GX,
                  double GY,
                  double alpha,
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G
                  )
{
    
    // calcualting F
    //boundary values
    for (int j=1; j <=jmax; j++)
    {
    F[0][j]=U[0][j];
    F[imax][j]=U[imax][j];
    }

    
// grid points
    int i,j=1;
    for(i = 1; i <= imax-1; i++) {
        for(j = 1; j<=jmax; j++) {
            F[i][j]=U[i][j]+ dt*(
            (1/Re)*(d2u_dx2( U[i-1][j], U[i][j], U[i+1][j],dx) + d2u_dx2( U[i][j-1], U[i][j], U[i][j+1], dy) )
            - du2_dx( U[i-1][j], U[i][j], U[i+1][j], dx, alpha )
            -duv_dx(V[i][j-1], V[i+1][j-1], V[i][j], V[i+1][j], U[i][j-1], U[i][j], U[i][j+1], dy, alpha)
            +GX );
        }
    }
    
    
    //calculating G
    
    
      //boundary values
    for (i=1; i <=imax; i++)
    {
    G[i][0]=V[i][0];
    G[i][jmax]=V[i][jmax];
    }

  
  // grid points  
    for(i = 1; i <= imax; i++) {
        for(j = 1; j<=jmax-1; j++) {
            G[i][j]=V[i][j]+ dt*(
                                 (1/Re)*(d2u_dx2( V[i-1][j], V[i][j], V[i+1][j], dx) + d2u_dx2( V[i][j-1], V[i][j], V[i][j+1], dy) )
            - du2_dx( V[i][j-1], V[i][j], V[i][j+1], dy, alpha )
            -duv_dx(U[i-1][j], U[i-1][j+1], U[i][j], U[i][j+1], V[i-1][j], V[i][j], V[i+1][j], dx, alpha)
            +GY );
        }
    }
    
    
}




//calculate rs
void calculate_rs(
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **F,
                  double **G,
                  double **RS
                  )
{
    int i,j=1;
    for(i = 1; i <= imax; i++) {
        for(j = 1; j<=jmax; j++) {
    RS[i][j] = (1/dt)*( du_dx(F[i-1][j],F[i][j],dx) + du_dx(G[i][j-1],G[i][j],dy) );
        }
    }
    }


        
//calculate uv
// can reduce looping effort by calculating for imax-1, jmax-1 once and then specifying end points in a different loop. not done right now.

        void calculate_uv(
                          double dt,
                          double dx,
                          double dy,
                          int imax,
                          int jmax,
                          double **U,
                          double **V,
                          double **F,
                          double **G,
                          double **P
                          )
        {
            //calcualting U
            
            int i,j=1;
            for(i = 1; i <= imax-1; i++) {
                for(j = 1; j<=jmax; j++) {
            U[i][j]= F[i][j]- dt*(du_dx(P[i][j],P[i+1][j],dx));
        }
            }
        
// calculating V
         for(i = 1; i <= imax; i++) {
                for(j = 1; j<=jmax-1; j++) {
   V[i][j]=G[i][j] - dt*(du_dx(P[i][j],P[i][j+1],dy));
      }
            }
            
        }



//calculate dt
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
  if(tau>0)
  {
    double coeff = Re/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  
    double min_value=fmin(coeff,fmin((dx/max_val_U(imax,jmax,U)),(dx/max_val_V(imax,jmax,V))));

    *dt=tau*min_value;
  }

}



	



// calculating discretised values
double du_dx( double uim1j, double uij, double dx) {
    double du_dx=(uij-uim1j)/dx;
    return du_dx;
}
// for pressure, the order is different (It is p(i+1,j)- p(ij) )

double d2u_dx2( double uim1j, double uij, double uip1j, double dx) {
    double d2u_dx2=(uim1j+uip1j-2*uij)/(dx*dx);
    return d2u_dx2;
}

double duv_dx (double uim1j, double uim1jp1, double uij, double uijp1, double vim1j, double vij, double vip1j, double dx, double alpha) {
    double duv_dx=( 1/(4*dx) )*( (uij+uijp1)*(vij+vip1j)- (uim1j+uim1jp1)*(vim1j+vij) +
                                alpha*( (fabs(uij+uijp1))*(vij-vip1j) - (fabs(uim1j+uim1jp1))*(vim1j-vij) ) );
    return duv_dx;
}






double du2_dx (double uim1j, double uij, double uip1j, double dx, double alpha) {
    double du2_dx= (1/(4*dx))*( (uij+uip1j)*(uij+uip1j) - (uim1j+uij)*(uim1j+uij) +
                               alpha*( (fabs(uij +uip1j))*(uij-uip1j) - (fabs(uim1j+ uij))*(uim1j-uij) ) );
    return du2_dx;
}


double max_val_U(int imax,int jmax,double** U)
{
  double max_val=-INFINITY;

  for (int i=0;i<=imax;i++)
  {
    for (int j=0;j<=jmax+1;j++)
       {
          if (max_val<fabs(U[i][j])) 
          max_val=fabs(U[i][j]);
       }
  }
  return max_val;
}

double max_val_V(int imax,int jmax,double** V)
{
  double max_val=-INFINITY;

  for (int i=0;i<=imax+1;i++)
  {
    for (int j=0;j<=jmax;j++)
       {
          if (max_val<fabs(V[i][j])) 
          max_val=fabs(V[i][j]);
       }
  }
  return max_val;
}

