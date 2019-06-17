#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include "boundary_val.h"
#include <stdio.h>
#include <time.h>

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

int main(int argn, char** args)
{
  clock_t start,end;
  double cpu_time_used;
  start=clock();
  
  double t=0;     
  int n=0;        /*t=time and n=iteration counter*/
  
  double Re,UI,VI,PI,GX,GY,t_end,
         xlength,ylength,dt,dx,dy,alpha,omg,tau,eps,dt_value;
  
  int imax,jmax,itermax; 
  
  read_parameters("cavity100.dat",&Re,&UI,&VI,&PI,&GX,&GY,&t_end,&xlength,&ylength,&dt,
                      &dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&dt_value);  
  
  /*Matrices of U,V,P are initialised with initial values*/
  double **U=NULL,**V=NULL,**P=NULL,**F=NULL,**G=NULL,**RS=NULL; 
  	
  /*Allocation of memory space to matrices*/ 
  P= matrix ( 0 , imax+1 , 0 , jmax+1 );
  U = matrix ( 0 , imax , 0 , jmax+1 );
  V = matrix ( 0 , imax+1 , 0 , jmax );
  F = matrix ( 0 , imax , 0 , jmax );
  G = matrix ( 0 , imax , 0 , jmax );
  RS= matrix ( 0 , imax+1 , 0 , jmax+1 );

  /*Creation of matrices with initial values*/
  init_uvp(UI,VI,PI,imax,jmax,U,V,P); 

  /*Initial Visualisation of U,V,P */

while (t<t_end)
{
  /*calculate_dt function*/
  calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);
   
  /*The boundary values are set.*/
  boundaryvalues(imax,jmax,U,V);

  /*Calculation of Fn and Gn*/
  calculate_fg(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,U,V,F,G);

  /*Calculation of RS*/
  calculate_rs(dt,dx,dy,imax,jmax,F,G,RS);

  int it=0;
  double d=(2*eps),*res=NULL;
  res=&d;

  while((it<itermax) && (*res>eps))
{
  sor(omg,dx,dy,imax,jmax,P,RS,res);
  it=it+1;
}
  
  /*Calculation of U and V for n+1 step*/
  calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P);

  /*Updation*/
  t=t+dt;
  n=n+1;

}
  /*Visualisation of U,V,P */
  write_vtkFile("results",n,xlength,ylength,imax,jmax,dx,dy,U,V,P);
  
  /*Memory free*/   
  free_matrix(P, 0 , imax+1 , 0 , jmax+1 );
  free_matrix (U, 0 , imax , 0 , jmax+1 );
  free_matrix (V, 0 , imax+1 , 0 , jmax );
  free_matrix (F, 0 , imax , 0 , jmax );
  free_matrix (G, 0 , imax , 0 , jmax );
  free_matrix (RS, 0 , imax+1 , 0 , jmax+1 );
   
  end=clock();
  cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
  printf("\n Total elapsed time: %f\n",cpu_time_used);

  return -1;
}
