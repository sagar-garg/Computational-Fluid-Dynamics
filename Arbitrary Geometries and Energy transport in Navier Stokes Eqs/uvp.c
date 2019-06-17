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
#include "init.h"
#include "uvp.h"



// declaring functions used for discretization
double du_dx( double uim1j, double uij, double dx);
double d2u_dx2( double uim1j, double uij, double uip1j, double dx);

double duv_dx (double uim1j, double uim1jp1, double uij, double uijp1, double vim1j, double vij, double vip1j, double dx, double alpha);
double du2_dx (double uim1j, double uij, double uip1j, double dx, double alpha);

double duT_dx (double uim1j, double uij, double Tim1j, double Tij, double Tip1j, double dx, double alpha);


//helper function for calculate_dt. find maximum value in a matrix
double max_val_U(int imax,int jmax,double** U);
double max_val_V(int imax,int jmax,double** V);


void calculate_fg(
                  double Re,
                  double GX,
                  double GY,
                  double alpha,
                  double beta,
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G,
                  double **T,
                  int **flag
                  )
{
    
    // calcualting F
    //boundary values
    for (int j=1; j <=jmax; j++)
    {
    F[0][j]=U[0][j];
    F[imax][j]=U[imax][j];
    
    }
    for(int i = 1; i < imax; i++){
        for (int j = 1; j <= jmax; j++)
        {
            switch(flag[i][j]){
                case 130: F[i-1][j] = U[i-1][j];
                    break;
                case 258: F[i][j] = U[i][j];
                    break;
                case 290: F[i][j] = U[i][j];
                    break;
                case 162: F[i-1][j] = U[i-1][j];
                    break;
                case 194: F[i-1][j] = U[i-1][j];
                    break;
                case 322: F[i][j] = U[i][j];
                    break;
                default:
                    break;
                    /*
                    
                     BW 130
                     BO 258
                     
                     BNO 290
                     BNW 162
                     BWS 194
                     BSO 322
                     */
            }
        }
    }

    //int imax_count = 0;
// grid points
    int i,j=1;
    for(i = 1; i <= imax-1; i++) {
        for(j = 1; j<=jmax; j++) {
            if(flag[i][j]==1){
            
            //imax_count++;
            
            F[i][j]=U[i][j]+ dt*(
            (1/Re)*(d2u_dx2( U[i-1][j], U[i][j], U[i+1][j],dx) + d2u_dx2( U[i][j-1], U[i][j], U[i][j+1], dy) )
            - du2_dx( U[i-1][j], U[i][j], U[i+1][j], dx, alpha )
            -duv_dx(V[i][j-1], V[i+1][j-1], V[i][j], V[i+1][j], U[i][j-1], U[i][j], U[i][j+1], dy, alpha)
            +GX )-(((beta*dt*GX)/2)*(T[i][j]+T[i+1][j]));

            //F[i][j]=F[i][j]-(((beta*dt*GX)/2)*(T[i][j]+T[i+1][j])); //Temperature dependent problem
            }
        }
    }

    
    
    
    //calculating G
    
    
      //boundary values
    for (i=1; i <=imax; i++)
    {
    G[i][0]=V[i][0];
    G[i][jmax]=V[i][jmax];
    
    }
        
    for(int i = 1; i <= imax; i++){
        for (int j = 1; j < jmax; j++)
            {
                switch(flag[i][j]){
                    case 34: G[i][j] = V[i][j];
                        break;
                    case 66: G[i][j-1] = V[i][j-1];
                        break;
                    case 290: G[i][j] = V[i][j];
                        break;
                    case 162: G[i][j] = V[i][j];
                        break;
                    case 194: G[i][j-1] = V[i][j-1];
                        break;
                    case 322: G[i][j-1] = V[i][j-1];
                        break;
                    default:
                        break;
                        
                    /*
                         BN 34
                         BS 66

                         BNO 290
                         BNW 162
                         BWS 194
                         BSO 322
                    */
                }
        }
    }

    //int jmax_count = 0;
  // grid points  
    for(i = 1; i <= imax; i++) {
        for(j = 1; j<=jmax-1; j++) {
            if(flag[i][j]==1){
                
            //jmax_count++;
            
            G[i][j]=V[i][j]+ dt*((1/Re)*(d2u_dx2( V[i-1][j], V[i][j], V[i+1][j], dx) + d2u_dx2( V[i][j-1], V[i][j], V[i][j+1], dy) )
            - du2_dx( V[i][j-1], V[i][j], V[i][j+1], dy, alpha )
            -duv_dx(U[i-1][j], U[i-1][j+1], U[i][j], U[i][j+1], V[i-1][j], V[i][j], V[i+1][j], dx, alpha)
            +GY)-(((beta*dt*GY)/2)*(T[i][j]+T[i][j+1]));

            //G[i][j]=G[i][j]-(((beta*dt*GY)/2)*(T[i][j]+T[i][j+1])); //Temperature dependent problem

            }
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
                  double **RS,
                  int **flag
                )
{
    int i,j=1;
    for(i = 1; i <= imax; i++) {
        for(j = 1; j<=jmax; j++) {
          if(flag[i][j]==1){
    RS[i][j] = (1/dt)*( du_dx(F[i-1][j],F[i][j],dx) + du_dx(G[i][j-1],G[i][j],dy) );
          }
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
                          double **P,
                          int **flag
                          )
        {
            //calcualting U
            
            int i,j=1;
            for(i = 1; i <= imax-1; i++) {
                for(j = 1; j<=jmax; j++) {
                    if(flag[i][j]==1){
            U[i][j]= F[i][j]- dt*(du_dx(P[i][j],P[i+1][j],dx));
        }
            }
            }
        
// calculating V
         for(i = 1; i <= imax; i++) {
                for(j = 1; j<=jmax-1; j++) {
                    if(flag[i][j]==1){
   V[i][j]=G[i][j] - dt*(du_dx(P[i][j],P[i][j+1],dy));
      }
            }
            
        }
        }



//calculate dt
void calculate_dt(
    const char *szFileName,
  double Re,
  double Pr,
  double *dt,
  double tau,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
)
{
//double dt=0;
if(tau<0.0)
    {
        READ_DOUBLE( szFileName, *dt  );
    }
else
    {
        double c= (1.0/pow(dx,2.0))+(1.0/pow(dy,2.0));
        double umax=max_ab(U,imax,jmax);
        double vmax=max_ab(V,imax,jmax);
        double min_1,min_2,min_3,min_last;
        // new stability constraint 
        min_3=((Re*Pr)/2.0)*(pow(c,-1));

        min_1=min_f(Re/(2.0*c),dx/umax);
        min_2=min_f(min_1,dy/vmax);
        min_last=min_f(min_2,min_3);
        *dt= tau*min_last;
    }
}
 /*   
double dt=0;

  if(tau>0)
  {
    double coeff = (Pr*Re)/(2.0  *  (  (1.0/(dx*dx))  +  (1.0/(dy*dy))   ));
  
    double min_value=fmin(coeff,fmin((dx/max_val_U(imax,jmax,U)),(dy/max_val_V(imax,jmax,V))));

    dt=tau*min_value;
  }
  return dt;

}
*/
/*
//Function which solves Energy equation to get Temperature value at n+1 time step
void calculate_temp(double dt,
                          double dx,
                          double dy,
                          int imax,
                          int jmax, double alpha, double Re, double Pr, double **T, double **U,
                          double **V)
{
    int i,j=1;//Calculation of inner grid points (Boindaries are already set)
    for(i = 1; i <= imax; i++) {
        for(j = 1; j<=jmax; j++) {
    T[i][j]=T[i][j]+ dt*((1/(Re*Pr))*(d2u_dx2( T[i-1][j], T[i][j], T[i+1][j], dx) + d2u_dx2( T[i][j-1], T[i][j], T[i][j+1], dy) )
    -(duv_dx(U[i-1][j], U[i-1][j], U[i][j], U[i][j], T[i-1][j], T[i][j], T[i+1][j], dx, alpha))
    -(duv_dx(V[i][j-1], V[i][j-1], V[i][j], V[i][j], T[i][j-1], T[i][j], T[i][j+1], dy, alpha)));
        }
    }

    //U[i][j+1]=U[i][j]=T[i][j]
    //U[i-1][j+1]=U[i-1][j]=T[i-1][j]
}

*/

//Function which solves Energy equation to get Temperature value at n+1 time step
void calculate_temp(double dt,
                          double dx,
                          double dy,
                          int imax,
                          int jmax, double alpha, double Re, double Pr, double **T, double **U,
                          double **V)
                          {
double **dumb=NULL; 
//init_matrix( dumb, 0, imax+1, 0 , jmax+1, 0 );
dumb = matrix ( 0 , imax+1 , 0 , jmax+1 );
  
for (int i=1; i<=imax; i++)
{ //calculating for inner grid points. boundary already set
    
    for (int j=1; j<=jmax; j++)
    {
    dumb[i][j]=T[i][j]+ 
    dt*( (1/(Re*Pr))*(d2u_dx2( T[i-1][j], T[i][j], T[i+1][j], dx) + d2u_dx2( T[i][j-1], T[i][j], T[i][j+1], dy) )
   -duT_dx(U[i-1][j], U[i][j], T[i-1][j], T[i][j], T[i+1][j], dx, alpha)
   -duT_dx(V[i][j-1], V[i][j], T[i][j-1], T[i][j], T[i][j+1], dy, alpha ));
    }
}

//updates have to happen after all values done.
for (int i=1; i<=imax; i++)
{
    for (int j=1; j<=jmax; j++)
    {
        //if cell=fluid cell
        T[i][j]=dumb[i][j];   
    }
}
 free_matrix (dumb, 0 , imax+1 , 0 , jmax+1 );
}



// calculating discretised values
double duT_dx (double uim1j, double uij, double Tim1j, double Tij, double Tip1j, double dx, double alpha) {
double duT_dx=(1/(2*dx))*( uij*(Tij+Tip1j) -uim1j*(Tim1j+Tij) + alpha*( (fabs(uij))*(Tij-Tip1j) - (fabs(uim1j))*(Tim1j-Tij) ) );
return duT_dx;
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

/*
double max_val_U(int imax,int jmax,double** U)
{
  double max_val=-INFINITY;

  for (int i=1;i<=(imax-1);i++)
  {
    for (int j=1;j<=jmax;j++)
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

  for (int i=1;i<=imax;i++)
  {
    for (int j=1;j<=(jmax-1);j++)
       {
          if (max_val<fabs(V[i][j])) 
          max_val=fabs(V[i][j]);
       }
  }
  return max_val;
}

*/