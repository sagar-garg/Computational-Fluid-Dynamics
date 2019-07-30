#include "helper.h"
#include "init.h"

int read_parameters(const char *szFileName,       /* name of the file */
                    int *num_boundary_particles,  /* number of particles present on the boundary */
                    int *num_inner_particles,     /* number of particles present in the interial */
                    double *rho,                  /* density of the particles */
                    double *radius,            /* radius of influence domain */ 
                    double *k_scale,            /* scaling factor depending on the kernel */
                    double *x_length,           /* length of the domain x-dir.*/
                    double *y_length,           /* length of the domain y-dir.*/
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                 /* initial pressure */ 
                    double *U_in,              /* U_velocity inflow */
                    double *V_in,              /* V_velocity inflow */
                    int *temp,                 /* Temp flag set to one for temperature dependent problem */
                    double *TI,                /* Initial temperature value */
                    double *Te,                /* East side temperature */
                    double *Tw,                /* West side temperature */
                    double *Tn,                /* North side temperature */
                    double *Ts,                /* South side temperature */   
                    double *dt,                /* time step */
                    int  *itermax,             /* max. number of iterations  */
                    double *eps,               /* accuracy bound for pressure*/
		                double *dt_value,          /* time for output */
                    int *iter_num,              /* printing after iteration */
                    double *gamma,             /* constant */
                    double *rho_reference)     /* density at zero time */
{
  READ_INT   ( szFileName, *num_boundary_particles );
  READ_INT   ( szFileName, *num_inner_particles );
  READ_DOUBLE( szFileName, *rho );
  READ_DOUBLE( szFileName, *radius );
  READ_DOUBLE( szFileName, *k_scale );
  READ_DOUBLE( szFileName, *x_length );
  READ_DOUBLE( szFileName, *y_length );
  READ_DOUBLE( szFileName, *Re );
  READ_DOUBLE( szFileName, *UI );
  READ_DOUBLE( szFileName, *VI );
  READ_DOUBLE( szFileName, *PI );
  READ_DOUBLE( szFileName, *U_in );
  READ_DOUBLE( szFileName, *V_in );
  READ_INT   ( szFileName, *temp );
  READ_DOUBLE( szFileName, *TI );
  READ_DOUBLE( szFileName, *Te );
  READ_DOUBLE( szFileName, *Tw );
  READ_DOUBLE( szFileName, *Tn );
  READ_DOUBLE( szFileName, *Ts );
  READ_DOUBLE( szFileName, *dt );    
  READ_INT   ( szFileName, *itermax );
  READ_DOUBLE( szFileName, *eps ); 
  READ_DOUBLE( szFileName, *dt_value );
  READ_INT   ( szFileName, *iter_num );
  READ_DOUBLE( szFileName, *gamma );
  READ_DOUBLE( szFileName, *rho_reference );

  return 1;
}

void init(
  struct particle* domainparticles,
  int num_particles,
  double* old_temperature, 
  double* old_u_velocity, 
  double* old_v_velocity,
  double UI,
  double VI,
  double TI)
{
  for (int i=0; i< num_particles; i++)
  {
    domainparticles[i].neighbour_distance = (double*) calloc(num_particles,sizeof(double)); //initialistion of neighbour_distance array
    domainparticles[i].weights = (double*) calloc(num_particles,sizeof(double)); //initialistion of neighbour_distance array
    domainparticles[i].dwdq = (double*) calloc(num_particles,sizeof(double)); //initialistion of dwdx array
    domainparticles[i].dwdx = (double*) calloc(num_particles,sizeof(double)); //initialistion of dwdx array
    domainparticles[i].dwdy = (double*) calloc(num_particles,sizeof(double)); //initialistion of dwdx array
    
    old_temperature[i] = TI;
    old_u_velocity[i] = UI;
    old_v_velocity[i] = VI;
    
  }




}


void set_input_inner(struct particle* domainparticles, double x_len, double y_len, int num_boundary_particles, int num_inner_particles,
double TI, double UI, double VI, double rho)
{   
int m  = num_inner_particles; // number of (inner) particles is perfect square
  int div = sqrt(m);  // if n = 81, div = 9, x_len = 10. 
  //Shouldnt it be n=100, div=10, x_len=10??
  double hx = x_len/div; // spacing in x-direction
  double hy = y_len/div; // spacing in y-direction
  int count = 0;
  for(int j = 0; j < div; j++){
    for(int i = 0; i< div; i++){
        domainparticles[count].temperature = TI; //
        domainparticles[count].rho = rho; // 
        domainparticles[count].mass = rho*hx*hy; 
      domainparticles[count].u_vel = UI ;
      domainparticles[count].v_vel = VI ;
      domainparticles[count].dvxdt = 0;
      domainparticles[count].dvydt = 0;
      domainparticles[count].drhodt = 0;

        
        count = count + 1;
  }
  }
}


  
  //count=m
 // // we go anti-clockwise: south east north west
void set_input_boundary(struct particle* domainparticles, double x_len, double y_len, int num_boundary_particles, int num_inner_particles,
double TE, double TW, double TN, double TS, double rho, double U_in, double V_in)
{   
  int m  = num_inner_particles; // number of (inner) particles is perfect square
  int div = sqrt(m);  // if n = 81, div = 9, x_len = 10. 
  //Shouldnt it be n=100, div=10, x_len=10??
  double hx = x_len/div; // spacing in x-direction
  double hy = y_len/div; // spacing in y-direction
  
  //printf("%f", hx);
  int count = m; //NOTE!!
  
  
  
  for(int j = 0; j < div; j++){ //south
      domainparticles[count].temperature = TS; //
      domainparticles[count].rho = rho; // 
      domainparticles[count].mass = rho*hx*hy; 
      domainparticles[count].u_vel = 0 ;
      domainparticles[count].v_vel = 0 ;
      domainparticles[count].dvxdt = 0;
      domainparticles[count].dvydt = 0;
      domainparticles[count].drhodt = 0;
      count = count + 1;
  }
  

  for(int j = 0; j < div; j++){ //east
      domainparticles[count].temperature = TE; //
      domainparticles[count].rho = rho; // 
      domainparticles[count].mass = rho*hx*hy; 
      domainparticles[count].u_vel = 0 ;
      domainparticles[count].v_vel = 0 ;
      domainparticles[count].dvxdt = 0;
      domainparticles[count].dvydt = 0;
      domainparticles[count].drhodt = 0;
      count = count + 1;
  }

  for(int j = 0; j < div; j++){ //north
      domainparticles[count].temperature = TN; //
      domainparticles[count].rho = rho; // 
      domainparticles[count].mass = rho*hx*hy; 
      domainparticles[count].u_vel = U_in ;
      domainparticles[count].v_vel = V_in ;
      domainparticles[count].dvxdt = 0;
      domainparticles[count].dvydt = 0;
      domainparticles[count].drhodt = 0;
      count = count + 1;
  }

  for(int j = 0; j < div; j++){ //west
      domainparticles[count].temperature = TW; //
      domainparticles[count].rho = rho; // 
      domainparticles[count].mass = rho*hx*hy; 
      domainparticles[count].u_vel = 0 ;
      domainparticles[count].v_vel = 0 ;
      domainparticles[count].dvxdt = 0;
      domainparticles[count].dvydt = 0;
      domainparticles[count].drhodt = 0;
      count = count + 1;
  }
  
}


