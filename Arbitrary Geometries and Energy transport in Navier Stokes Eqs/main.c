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

  int number=0;   //For problem number
  int temp=0;     //Temp Flag, ON if problem is temperature dependent
  const char *problem=NULL;  //Problem File name
  printf("Provide the problem name: Press\n 1 - The Karman Vortex street\n 2 - Flow over a step \n 3 - Natural Convection 1 \n 4 - Natural Convection 2 \n 5 - Fluid trap \n 6 - Rayleigh-Benard Convection\n 7 - Any other problem\n");
  scanf("%d", &number);
 
//For new problem
  char new_problem[100];
  int prob_num;
 
  switch (number)
  {
    case 1:
    problem="the-Karman-Vortex-street.dat";
    temp=0;
    break;
    
    case 2:
    problem="flow-over-a-step.dat";
    temp=0;
    break;
    
    case 3:
    problem="natural-convection-1.dat";
    temp=1;
    break;
    
    case 4 :
    problem="natural-convection-2.dat";
    temp=1;
    break;
    
    case 5: 
    problem="fluid-trap.dat";
    temp=1;
    break;
    
    case 6:
    problem="rayleigh-benard-convection.dat";
    temp=1;
    break;
    
    case 7://problem=Any other problem;
    {
    printf("Is it Temperature dependent or independent problem?\n Press 1 : Temperature dependent \n Press 0 : Temperature independent\n");
    scanf("%d",&prob_num);
    switch (prob_num)
    {
      case 1:
      temp=1;//Temperature dependent
      break;
      case 0:
      temp=0;//Temperature independent
      break;
      default:
      printf("\nError: Enter 1 or 0");
      break;
    }

    printf("\nEnter the name of dat file (with .dat extension):");
    scanf("%s",new_problem);
    problem=new_problem;
    break;
    }

    default:
    printf("Error in number pressed");
    break;
  }
  
  
  char geometry[100];//Geometry File name (need to read from Dat file)Just to initialise
  int img_cnt=1;    //For counting images of VTK file
  double t=0;      //t=time and
  int n=0;         //n=iteration counter

  double Re,Pr,UI,VI,U_in,V_in,PI,TI,Te,Tw,Tn,Ts,GX,GY,t_end,
         xlength,ylength,dt,dx,dy,alpha,omg,tau,beta,eps,dt_value;
  
  int imax,jmax,itermax; 

  read_parameters(problem,geometry,&Re,&Pr,&UI,&VI,&U_in,&V_in,&PI,&TI,&Te,&Tw,&Tn,&Ts,&GX,&GY,&t_end,&xlength,&ylength,&dt,
                      &dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&beta,&itermax,&eps,&dt_value);  

  int** flag=NULL;
  flag = imatrix(0, (imax+2)-1, 0, (jmax+2)-1); //Initiation of  Flag matrix.
  int count=init_flag(problem,geometry,imax,jmax,flag); //Setting flags to each cell.//Count is the number of fluid cells 
  
  //printing Flag matrix
  /*
  for(int i=0;i<imax+2;i++)
  {
    printf("\n");
    for(int j=0;j<jmax+2;j++)
    {
      printf("%d ",flag[i][j]);
    }
  }
  printf("\n");
  */
  
  /*Matrices of U,V,P are initialised with initial values*/
  double **U=NULL,**V=NULL,**P=NULL,**T=NULL,**F=NULL,**G=NULL,**RS=NULL; 
  	
  /*Allocation of memory space to matrices*/ 
  P = matrix ( 0 , imax+1 , 0 , jmax+1 );
  U = matrix ( 0 , imax , 0 , jmax+1 );
  V = matrix ( 0 , imax+1 , 0 , jmax );
  T = matrix ( 0 , imax+1 , 0 , jmax+1 );
  F = matrix ( 0 , imax , 0 , jmax );
  G = matrix ( 0 , imax , 0 , jmax );
  RS = matrix ( 0 , imax+1 , 0 , jmax+1 );

  /*Creation of matrices with initial values*/
  init_uvp(UI,VI,PI,TI,imax,jmax,U,V,P,T);
 
  while (t<t_end)
{
  /*calculate_dt function Updated one addition of Pr.*/ 
  calculate_dt(problem,Re,Pr,&dt,tau,dx,dy,imax,jmax,U,V);

  /*The boundary values are set.*/ 
  if(temp!=1){
    //Temperature independent Problem 
    boundaryvalues(imax,jmax,U,V,flag,UI,VI,U_in,V_in);
    spec_boundary_val(imax,jmax,U,V,flag,UI,VI);
  }
  else{ 
    //Temperature dependent problem
    boundaryvalues_t(imax,jmax,U,V,T,flag,UI,VI,U_in,V_in,Te,Tw,Tn,Ts);
    spec_boundary_val_t(imax,jmax,U,V,flag,T);
    
    //Compute T(n+1) using T(n) using expl. Euler and Energy equation
    calculate_temp(dt,dx,dy,imax,jmax,alpha,Re,Pr,T,U,V);
  } 

  //Visualisation at t=0, Updated for Temperature also , But what happen for temperature independent problem 
  if(n==0)
  {
    write_vtkFile(problem,n,xlength,ylength,imax,jmax,temp,dx,dy,U,V,P,T);//Temperature added
  }

  /*Calculation of Fn and Gn*/ 
  //updated according to new formulae of F and G (Temperature consideration)
  calculate_fg(Re,GX,GY,alpha,beta,dt,dx,dy,imax,jmax,U,V,F,G,T,flag);

  /*Calculation of RS*/
  calculate_rs(dt,dx,dy,imax,jmax,F,G,RS,flag);

  int it=0;
  double d=(2*eps),*res=NULL;
  res=&d;

  while((it<itermax) && (*res>eps))
{
  sor(omg,dx,dy,imax,jmax,P,RS,res,flag,count);//Count is the number of fluid cells 
  it=it+1;
}
 // printf("dt= %lf, loop= %d, SOR_iter= %d \n",t,n,it);
 printf("dt= %lf, SOR_iter= %d \n",t,it);

if (it==itermax){
  printf ("Warning: max iteration reached. Solution won't converge in this step.\n");
}

  //Already updated according to the temperature (As We have updated F AND G)
  /*Calculation of U and V for n+1 step*/
  calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,flag);
  
  /*Visualisation of U,V,P */
  int num_file = img_cnt*dt_value;
  if (t>=num_file)
{
  write_vtkFile(problem,n,xlength,ylength,imax,jmax,temp,dx,dy,U,V,P,T);//Temperature added
  img_cnt=img_cnt+1;
}
  t=t+dt;
  n=n+1;
}
 //write_vtkFile(problem,n,xlength,ylength,imax,jmax,temp,dx,dy,U,V,P,T);//Temperature added

  /*Memory free*/   
  free_matrix (P, 0 , imax+1 , 0 , jmax+1 );
  free_matrix (U, 0 , imax , 0 , jmax+1 );
  free_matrix (V, 0 , imax+1 , 0 , jmax );
  free_matrix (F, 0 , imax , 0 , jmax );
  free_matrix (G, 0 , imax , 0 , jmax );
  free_matrix (RS, 0 , imax+1 , 0 , jmax+1 );
  free_matrix (T, 0 , imax+1 , 0 , jmax+1 );
  
  end=clock();
  cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
  printf("\n Total elapsed time: %f\n",cpu_time_used);

  return -1;
}
