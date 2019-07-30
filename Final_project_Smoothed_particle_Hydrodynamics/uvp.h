#ifndef __UVP_H__
#define __UVP_H__

void calculate_temp_laplacian(struct particle* domainparticles, int num_particles, double h,double k_scale, double* temp_laplacian);

void calculate_temperature(struct particle* domainparticles, double* old_temperature, double* temp_laplacian, double dt, int num_particles);

void calculate_density(struct particle* domainparticles, int num_particles, double dt); 

void pressure_update(struct particle* domainparticles, int num_particles, double ref_density, double gamma);

void calculate_velocity(struct particle* domainparticles, int num_inner_particles, double dt, double h);

void calculate_position(struct particle* domainparticles, int num_inner_sparticles, double dt);

double calculate_dt(double h);

#endif

