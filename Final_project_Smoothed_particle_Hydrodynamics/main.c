#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "helper.h"
#include "init.h"
#include "uvp.h"
#include "visual.h"

int main(int argc, char* argv[]) 
{
  clock_t start,end;
  double cpu_time_used;
  start=clock();

  int number=0; //For problem number
  const char *problem=NULL;  //Problem File name
  printf("Provide the problem name: Press\n 1 - Heated_plate \n 2 - Lid driven cavity (visualisation of 9 particles) \n 3 - Lid driven cavity (visualisation of all particles)\n");
  scanf("%d", &number);

  switch (number)
  {
    case 1:
    problem="heated_plate.dat";
    break;
    
    case 2:
    problem="lid_driven_cavity_9_particles.dat";
    break;
    
    case 3:
    problem="lid_driven_cavity_all_particles.dat";
    break;

    default:
    printf("Error in number pressed");
    break;
  }

  char sol_directory[100];
  char sol_folder[100];
  struct stat st;
  sprintf(sol_folder,"solution_%s",problem);
  if(stat(sol_folder,&st)==-1)
    mkdir(sol_folder,0700);

  sprintf(sol_directory,"solution_%s/sol",problem);

  int num_boundary_particles, num_inner_particles, temp, itermax, iter_num;
  double rho, radius, k_scale, x_length, y_length, Re, UI, VI, PI, U_in, V_in, TI, Te, Tw, Tn, Ts, dt, eps, dt_value, gamma, rho_reference ;

  //Reading all parameters from DAT file
  read_parameters(problem,&num_boundary_particles,&num_inner_particles,&rho,&radius,&k_scale,&x_length,&y_length,&Re,&UI,&VI,&PI,&U_in,&V_in,&temp,&TI,&Te,&Tw,&Tn,&Ts,&dt,&itermax,&eps,&dt_value,&iter_num,&gamma,&rho_reference);

  int num_particles = num_boundary_particles+num_inner_particles; //total number of particles
  double h = radius/k_scale;
  //printf("%d %d, %d", num_boundary_particles, num_inner_particles, num_particles);
  
  struct particle* domainparticles; 
  domainparticles= (struct particle*) calloc(num_particles,sizeof(struct particle));
  
  double* old_temperature; //temprory array to store the temerature of each particle
  old_temperature= (double*) calloc (num_particles,sizeof(double));

  double* temp_laplacian; //array to store the laplacian of temerature of each particle
  temp_laplacian= (double*) calloc (num_particles,sizeof(double));

  double* old_u_velocity; //array to store the laplacian of temerature of each particle
  old_u_velocity = (double*) calloc (num_particles,sizeof(double));

  double* old_v_velocity; //array to store the laplacian of temerature of each particle
  old_v_velocity= (double*) calloc (num_particles,sizeof(double));

  //Initialises x and y coordinates of all particles 
  set_coordinates(domainparticles, x_length, y_length, num_boundary_particles,num_inner_particles);

  //intialisation of neighbour_distance, weights, dwdx, dwdy arrays 
  init(domainparticles,num_particles,old_temperature, old_u_velocity, old_v_velocity, UI, VI, TI);

  //initialize inputs of T,P, rho, m, etc. on inner cells and boundary
  set_input_inner(domainparticles, x_length, y_length, num_boundary_particles, num_inner_particles, TI, UI, VI, rho);

  printf("\n %f  %f \n ", domainparticles[2].rho, rho_reference);
  double max_error=2*eps; 
  int it_cnt=1;

  while (max_error >= eps  &&  it_cnt < itermax) //Accuracy achieved in current time step
  {
    
    
    for (int i=0; i<num_particles; i++)
    {
      for(int j=0; j<num_particles ; j++)
      {
        domainparticles[i].neighbour_distance[j]=calculate_distance(domainparticles[i].x_coordinate,domainparticles[i].y_coordinate,domainparticles[j].x_coordinate,domainparticles[j].y_coordinate);
      }
    }

    //Run the kernel for each particle 
    for (int i=0; i<=num_particles-1; i++)
    {
      double sum = 0;
      for(int j=0; j<=num_particles-1 ; j++)
      {
        domainparticles[i].weights[j]=kernel(domainparticles[i].neighbour_distance[j],h);
        sum += domainparticles[i].weights[j];
        domainparticles[i].dwdq[j]=kernel_grad(domainparticles[i].neighbour_distance[j],h);
        domainparticles[i].dwdx[j]=kernel_grad_x(domainparticles[i].neighbour_distance[j],h, domainparticles[i].x_coordinate, domainparticles[j].x_coordinate);
        domainparticles[i].dwdy[j]=kernel_grad_x(domainparticles[i].neighbour_distance[j],h, domainparticles[i].y_coordinate, domainparticles[j].y_coordinate);
      }
      for(int j=0; j<=num_particles-1 ; j++)
      {
        domainparticles[i].weights[j] = domainparticles[i].weights[j]/sum;
      }
    }

  
    if(temp==1) //Temperature dependent problem
    {
      // Boundary need to change according to the Bouncing prociple
      set_input_boundary(domainparticles, x_length, y_length, num_boundary_particles, num_inner_particles, Te, Tw, Tn, Ts, rho, U_in, V_in); 

      //Temperature Laplacian
      calculate_temp_laplacian(domainparticles, num_particles, h, k_scale,temp_laplacian);

      //initially storing the temperature values of all the particle inside the array 
      for(int i=0; i<num_particles; i++)
      {
        old_temperature[i]=domainparticles[i].temperature;
      }

      //Temporal integration  
      calculate_temperature(domainparticles, old_temperature, temp_laplacian, dt, num_inner_particles);

      //Check convergence 
      check_convergence(domainparticles, &max_error, num_inner_particles, num_boundary_particles, old_temperature, old_u_velocity, old_v_velocity,temp);
    }

    else // For cavity problem
    {
      // Boundary need to change according to the Bouncing prociple
      set_input_boundary(domainparticles, x_length, y_length, num_boundary_particles, num_inner_particles, Te, Tw, Tn, Ts, rho,U_in,V_in); 
      
      calculate_density(domainparticles, num_particles, dt);

      pressure_update(domainparticles, num_particles, rho_reference, gamma);
      
      //For Convegence
      for(int i=0; i<num_particles; i++)
      {
        old_u_velocity[i]=domainparticles[i].u_vel;
        old_v_velocity[i]=domainparticles[i].v_vel;
      }

      calculate_velocity(domainparticles, num_inner_particles, dt, h);

      calculate_position(domainparticles, num_inner_particles, dt);

      reflect(domainparticles, num_inner_particles, x_length, y_length, U_in);

      //Check convergence
      if(it_cnt > 10000)
      { 
        check_convergence(domainparticles, &max_error, num_inner_particles, num_boundary_particles, old_temperature, old_u_velocity, old_v_velocity,temp);
      }
    }
  
    it_cnt++;
    printf("\nIteration:%d - max_error:%lf",it_cnt,max_error);

    // Visualisation of U,V,P 
    if (it_cnt % iter_num == 0)
    {
      write_vtkFile(sol_directory,it_cnt,x_length,y_length,0,0,num_inner_particles,num_boundary_particles,(  x_length/sqrt(num_inner_particles)  ), (  y_length/sqrt(num_inner_particles)  ),domainparticles,temp); 
    }

  }

  if(it_cnt<itermax)
  {
      printf("\nsolution converged to given acccuracy at iteration = %d",it_cnt);
  }
  else
  {
      printf("\nsolution at ietmax = %d",it_cnt);
  }

  end=clock();
  cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
  printf("\nTotal elapsed time: %f\n",cpu_time_used);

  //Cleaning all the dynamic memory 
  cleaning(domainparticles,num_particles,old_temperature,old_u_velocity,old_v_velocity,temp_laplacian);
  return 0;
}
