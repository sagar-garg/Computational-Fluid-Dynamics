#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include "boundary_val.h"
#include <stdio.h>
#include <time.h>
#include "precice/SolverInterfaceC.h"
#include "precice_adapter.h"
#include <omp.h>
#include<string.h>

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
 */

int main(int argn, char** args)
{
  clock_t start,end;
  double cpu_time_used;
  start=clock();
  int count_num_fluid=0;
  int num_coupling_cells=0;
  int implicit=0; //If we want to solve given Coupling problem by Implicit COuling it sets to 1

  int number=0;   //For problem number
  //int temp=0;     //Temp Flag, ON if problem is temperature dependent
  const char *problem=NULL;  //Problem File name
  printf("Provide the problem name: Press\n 1 - Forced convection over a heated plate \n 2 - Natural convection in cavity \n 3 - F1-heat-exchange \n 4 - F2-heat-exchange \n 5 - Any other problem\n");
  scanf("%d", &number);
 
//For new problem
  char new_problem[100];
  int prob_num;
  double Th = 0; //required for dimention conversion
  double Tc = 0;  ////required for dimention conversion
 
  switch (number)
  {
    case 1:
    problem="heated-plate.dat";
    //temp=1;
    implicit=0;
    Th = 10; //max_T(temperatureCoupled, num_coupling_cells);
    Tc = 0;  //min_T(temperatureCoupled, num_coupling_cells);
    break;
    
    case 2:
    problem="convection.dat";
    //temp=1;
    implicit=1;
    Th = 10; 
    Tc = 0;  
    break;
    
    case 3:
    problem="F1-heat-exchange.dat";
    //temp=1;
    implicit=0;
    Th = 5; 
    Tc =0;  
    break;
    
    case 4 :
    problem="F2-heat-exchange.dat";
    //temp=1;
    implicit=0;
    Th = 5; 
    Tc =0;  
    break;
    
    case 5://problem=Any other problem . Yet to add the conndition of implicit coupling
    {
    printf("Is it Implicit or Explicit Coupling problem?\n Press 1 : Implicit \n Press 0 : Explicit\n");
    scanf("%d",&prob_num);
    switch (prob_num)
    {
      case 1:
      implicit=1;//
      break;
      case 0:
      implicit=0;//
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
  char precice_config[100], participant_name[100], mesh_name[100], read_data_name[100], write_data_name[100]; //Reading strings for PreCICE
  int img_cnt=1;    //For counting images of VTK file
  double t=0,t_cp=0;      //t=time and t_cp backup time required for implicit precice coupling
  int n=0;         //n=iteration counter

  double Re,Pr,UI,VI,U_in,V_in,PI,TI,Te,Tw,Tn,Ts,GX,GY,t_end,
         xlength,ylength,x_origin,y_origin,dt,dx,dy,alpha,omg,tau,beta,eps,dt_value;
  
  int imax,jmax,itermax; 

  read_parameters(problem,geometry,precice_config,participant_name,mesh_name,read_data_name,write_data_name,&Re,&Pr,&UI,&VI,
                  &U_in,&V_in,&PI,&TI,&Te,&Tw,&Tn,&Ts,&GX,&GY,&t_end,&xlength,&ylength,&x_origin,&y_origin,&dt,
                  &dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&beta,&itermax,&eps,&dt_value);  

  int** flag=NULL;
  flag = imatrix(0, (imax+2)-1, 0, (jmax+2)-1); //Initiation of Flag matrix.
  init_flag(problem,geometry,imax,jmax,flag,&count_num_fluid,&num_coupling_cells); //Setting flags to each cell.//Count is the number of fluid cells 
  printf("number = %d\n, fluid= %d \n", num_coupling_cells, count_num_fluid );

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
  double **U=NULL,**V=NULL,**P=NULL,**T=NULL,**F=NULL,**G=NULL,**RS=NULL,**U_cp=NULL, **V_cp=NULL, **T_cp=NULL;   
  	
  /*Allocation of memory space to matrices*/ 
  P = matrix ( 0 , imax+1 , 0 , jmax+1 );
  U = matrix ( 0 , imax , 0 , jmax+1 );
  V = matrix ( 0 , imax+1 , 0 , jmax );
  T = matrix ( 0 , imax+1 , 0 , jmax+1 );
  F = matrix ( 0 , imax , 0 , jmax );
  G = matrix ( 0 , imax , 0 , jmax );
  RS = matrix ( 0 , imax+1 , 0 , jmax+1 );
  
  /*Allocation of memory space to backup matrices*/
  U_cp = matrix ( 0 , imax , 0 , jmax+1 );
  V_cp = matrix ( 0 , imax+1 , 0 , jmax );
  T_cp=  matrix( 0 , imax+1 , 0 , jmax+1);
  
  //printf("Debugging\n");

  /*Creation of matrices with initial values*/
  init_uvp(UI,VI,PI,TI,imax,jmax,U,V,P,T,U_cp,V_cp,T_cp);
  
  /* Initialize PreCICE -constructor and  configuration */ 
  precicec_createSolverInterface(participant_name,precice_config,0,1);
  int dim = precicec_getDimensions();
  printf("dimensions of the problem is %d \n",dim);

  //Define coupling mesh
  int meshID = precicec_getMeshID(mesh_name);
  printf("Mesh id = %d\n",meshID);
  int* vertexIDs= precice_set_interface_vertices(imax, jmax, dx, dy, x_origin, y_origin, num_coupling_cells, meshID, flag); //Get couling cell ids
  
  /*for (int k=0 ; k<num_coupling_cells; k++)
  {
    printf("%d ",vertexIDs[k]);
  }
  printf("\n");
  */

  //define dirichlet part of coupling writtten by this solver
  int temperatureID = precicec_getDataID(write_data_name,meshID);
  double* temperatureCoupled = (double*) malloc(sizeof(double) * num_coupling_cells);
  //precice_write_temperature(imax, jmax, num_coupling_cells, temperatureCoupled, vertexIDs,temperatureID, T, flag);
  
  //define Newmann part of coupling read by this solver
  int heatFluxID = precicec_getDataID(read_data_name,meshID);
  double* heatfluxCoupled = (double*) malloc(sizeof(double)* num_coupling_cells);
  //set_coupling_boundary(num_coupling_cells,imax, jmax, dx,dy,heatfluxCoupled, T, flag);

  //Call precicec_initialize which will calculte the maxmimum time step to be used 
  double precice_dt= precicec_initialize();
  
  //const std::string& cowic = precice::constants::actionWriteIterationCheckpoint();
  const char* coric = precicec_actionReadIterationCheckpoint();
  const char* cowic = precicec_actionWriteIterationCheckpoint();
  

  //Initialized data at coupling interface
  precice_write_temperature(imax, jmax, num_coupling_cells, temperatureCoupled, vertexIDs,temperatureID, T, flag,Th,Tc);
  precicec_initialize_data();
  precicec_readBlockScalarData(heatFluxID, num_coupling_cells, vertexIDs, heatfluxCoupled);

//while (t<t_end)
while(precicec_isCouplingOngoing())//Time loop
{
 
  if(implicit==1)
  {
  if(precicec_isActionRequired(cowic))
  {
      write_checkpoint(t,U,V,T,&t_cp,U_cp,V_cp,T_cp,imax,jmax);
		  //precicec.fulfilledAction(cowic);
      precicec_fulfilledAction(cowic);
	}
  }

  /*calculate_dt function Updated one addition of Pr and precice_dt.*/ 
  calculate_dt(problem,Re,Pr,&dt,tau,dx,dy,imax,jmax,U,V,precice_dt);
  //dt = 0.00005; //testing if constant dt works
  //printf("dt = %lf\n",dt);
 
  //printf("\n\nBarrier\n\n");

  /*The boundary values are set.*/ 
 // if(temp!=1){
    //Temperature independent Problem 
   // boundaryvalues(imax,jmax,U,V,flag,UI,VI,U_in,V_in);
   // spec_boundary_val(imax,jmax,U,V,flag,UI,VI);
//  }
  //else{ 
    //Temperature dependent problem
    boundaryvalues_t(imax,jmax,U,V,T,flag,UI,VI,U_in,V_in,Te,Tw,Tn,Ts);
    spec_boundary_val_t(imax,jmax,U,V,flag,T); 
/*   
        for(int i=0;i<imax+2;i++)
  {
      printf("%lf ",T[i][0]);
    
  }
  printf("\n\n");
*/ 
    // This heat flux has to be Dimensionless 
    set_coupling_boundary(num_coupling_cells,imax, jmax, dx,dy,heatfluxCoupled, T, flag);//overwrite the Temp. values at coupling cells
 /*
    for(int i=0;i<imax+2;i++)
  {
      printf("%lf ",T[i][0]);
    
  }
  printf("\n\n");
  */

    //Compute T(n+1) using T(n) using expl. Euler and Energy equation
    calculate_temp(dt,dx,dy,imax,jmax,alpha,Re,Pr,T,U,V);
  //} 

  /*Visualisation at t=0, Updated for Temperature also , But what happen for temperature independent problem 
  if(n==0)
  {
    //write_vtkFile(problem,n,xlength,ylength,imax,jmax,temp,dx,dy,U,V,P,T);//Temperature added
  }*/

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
  sor(omg,dx,dy,imax,jmax,P,RS,res,flag,count_num_fluid);//Count is the number of fluid cells 
  it=it+1;
}
 printf("dt= %lf, SOR_iter= %d \n",dt,it);

if (it==itermax){
  printf ("Warning: max iteration reached. Solution won't converge in this step.\n");
}

  //Already updated according to the temperature (As We have updated F AND G)
  /*Calculation of U and V for n+1 step*/
  calculate_uv(dt,dx,dy,imax,jmax,U,V,F,G,P,flag);
/*
  for(int i=0;i<(num_coupling_cells);i++)
  {
      printf("%lf ",temperatureCoupled[i]);
    
  }
  printf("\n");
*/

  precice_write_temperature(imax, jmax, num_coupling_cells, temperatureCoupled, vertexIDs,temperatureID, T, flag,Th,Tc);
  
  
  precice_dt=precicec_advance(dt);
  //printf("Precise_dt = %lf\n",precice_dt);
  
  if(implicit==1)
  {
  if(precicec_isActionRequired(coric))
  {
      restore_checkpoint(&t,U,V,T,t_cp,U_cp,V_cp,T_cp,imax,jmax);
		  precicec_fulfilledAction(coric);
	}
  }


  precicec_readBlockScalarData(heatFluxID,num_coupling_cells,vertexIDs,heatfluxCoupled);
  double deltat = Th - Tc;
  for(int i = 0; i < num_coupling_cells; i++){
      if(deltat != 0){
      heatfluxCoupled[i] = heatfluxCoupled[i]/(100*deltat);
      }
    }

  //Here we need to convert dimensioned Heatflux to dimentionless Heatflux

  /*Visualisation of U,V,P */
  int num_file = img_cnt*dt_value;
  if (t>=num_file)
{
  write_vtkFile(problem,n,xlength,ylength,x_origin,y_origin,imax,jmax,1,dx,dy,U,V,P,T);//Temperature added
  img_cnt=img_cnt+1;
}
  t=t+dt;
  n=n+1;
}
/*
  for(int i=0;i<imax+2;i++)
  {
    printf("\n");
    for(int j=0;j<jmax+2;j++)
    {
      printf("%lf ",T[i][j]);
    }
  }
  printf("\n");
  */

  precicec_finalize();

 //write_vtkFile(problem,n,xlength,ylength,imax,jmax,temp,dx,dy,U,V,P,T);//Temperature added

  /*Memory free*/   
  free_matrix (P, 0 , imax+1 , 0 , jmax+1 );
  free_matrix (U, 0 , imax , 0 , jmax+1 );
  free_matrix (V, 0 , imax+1 , 0 , jmax );
  free_matrix (F, 0 , imax , 0 , jmax );
  free_matrix (G, 0 , imax , 0 , jmax );
  free_matrix (RS, 0 , imax+1 , 0 , jmax+1 );
  free_matrix (T, 0 , imax+1 , 0 , jmax+1 );
  free_matrix (T_cp, 0 , imax+1 , 0 , jmax+1 );
  free_matrix (U_cp, 0 , imax , 0 , jmax+1 );
  free_matrix (V_cp, 0 , imax+1 , 0 , jmax );
  free(temperatureCoupled);
  free(heatfluxCoupled);

  end=clock();
  cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
  printf("\n Total elapsed time: %f\n",cpu_time_used);

  return -1;
}
