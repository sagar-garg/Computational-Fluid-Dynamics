#ifndef __INIT_H_
#define __INIT_H_

/**
 * The arrays neighbour_distance and weights are initialized 
 */

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
                    double *PI,                /* velocity y-direction */
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
                    double *rho_reference);     /* density at zero time */

void init(
  struct particle* domainparticles,
  int num_particles,
  double* old_temperature, 
  double* old_u_velocity, 
  double* old_v_velocity,
  double UI,
  double VI,
  double TI);

void set_input_inner(struct particle* domainparticles, double x_len, double y_len, int num_boundary_particles, int num_inner_particles, 
double TI, double UI, double VI, double rho);

void set_input_boundary(struct particle* domainparticles, double x_len, double y_len, int num_boundary_particles, int num_inner_particles,
double TE, double TW, double TN, double TS, double rho, double U_in, double V_in);




#endif
