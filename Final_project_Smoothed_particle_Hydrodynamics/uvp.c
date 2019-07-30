#include <stdio.h>
#include <math.h>
#include "helper.h"
#include "init.h"
#include "uvp.h"


double calculate_dt(double h){
    double dt = 0;
    dt = (0.5*0.125*h*h*1000)/8.9;
    return dt;
}

void calculate_temp_laplacian(struct particle* domainparticles, int num_particles, double h, double k_scale, double* temp_laplacian)
{
    double sum;
    //double k_scale=2; what to do with this

    for (int i=0; i<=num_particles-1; i++)
    {
        sum=0;
        for(int j=0; j<=num_particles-1 ; j++)
        {
            //Checking condition of inner particles and particles in the influence domain only
            if(domainparticles[i].particle_id >= 0    &&   domainparticles[i].neighbour_distance[j] <= k_scale*h  &&   i!=j) 
            {
                sum=sum +
                    ( (domainparticles[j].mass/domainparticles[j].rho)  *  
                      (domainparticles[i].dwdq[j])  *  
                      (  1  /  domainparticles[i].neighbour_distance[j])  *
                      (domainparticles[i].temperature - domainparticles[j].temperature)   ) ;

                     
            }
        }
        //printf("\nsum = %lf",sum);
        temp_laplacian[i]=2*sum;
    }
    
}


void calculate_temperature(struct particle* domainparticles, double* old_temperature, double* temp_laplacian, double dt, int num_inner_particles)
{
    //finding new temperature using old temperature values and laplacian
    for(int i=0; i<num_inner_particles; i++)
    {
         // Calculation for inner particles only 
            domainparticles[i].temperature= old_temperature[i] + temp_laplacian[i]*dt;  //Euler is used
    }
}


void calculate_density(struct particle* domainparticles, int num_particles, double dt) 
{
    double dvx, dvy;
    for(int i = 0; i < num_particles; i++){
        for (int j = 0; j < num_particles; j++)
        {
            if(domainparticles[i].weights[j] > 0.0 && i!=j){
                dvx = domainparticles[i].u_vel - domainparticles[j].u_vel;
                dvy = domainparticles[i].v_vel - domainparticles[j].v_vel;
                domainparticles[i].drhodt = domainparticles[i].drhodt + (domainparticles[j].mass)*(dvx*domainparticles[i].dwdx[j] + dvy*domainparticles[i].dwdy[j]);
            }
        }
    domainparticles[i].rho = domainparticles[i].rho + dt*(domainparticles[i].drhodt);
    }
}


void pressure_update(struct particle* domainparticles, int num_particles, double ref_density, double gamma){
  for(int i = 0; i < num_particles; i++ ){
    domainparticles[i].pressure = ((0.85*100000)*(pow((domainparticles[i].rho/ref_density), gamma) - 1.0));
  }
}


void calculate_velocity(struct particle* domainparticles, int num_inner_particles, double dt, double h)
{
    double dvx, dvy, dx, dy;
    for(int i = 0; i < num_inner_particles; i++){
        for (int j = 0; j < num_inner_particles; j++)
        {
            if(domainparticles[i].weights[j] > 0.0 && i!=j){ //just the term1 of dvxdt, dvydt
      //          domainparticles[i].dvxdt -= (domainparticles[j].mass)*(((domainparticles[i].pressure)/(pow(domainparticles[i].rho, 2.0))) + ((domainparticles[j].pressure)/(pow(domainparticles[j].rho, 2.0))))*domainparticles[i].dwdx[j];
       //         domainparticles[i].dvydt -= (domainparticles[j].mass)*(((domainparticles[i].pressure)/(pow(domainparticles[i].rho, 2.0))) + ((domainparticles[j].pressure)/(pow(domainparticles[j].rho, 2.0))))*domainparticles[i].dwdy[j];
                domainparticles[i].dvxdt += (domainparticles[j].mass)*(((domainparticles[i].pressure)/(pow(domainparticles[i].rho, 2.0))) + ((domainparticles[j].pressure)/(pow(domainparticles[j].rho, 2.0))))*domainparticles[i].dwdx[j];
                domainparticles[i].dvydt += (domainparticles[j].mass)*(((domainparticles[i].pressure)/(pow(domainparticles[i].rho, 2.0))) + ((domainparticles[j].pressure)/(pow(domainparticles[j].rho, 2.0))))*domainparticles[i].dwdy[j];
     
            }
        }
    }
    double sumx = 0, sumy = 0;
    for(int i = 0; i < num_inner_particles; i++){ //for term2 of dvxdt, dvydt
        for (int j = 0; j < num_inner_particles; j++)
        {
            if(domainparticles[i].weights[j] > 0.0 && i!=j){
                dvx = domainparticles[i].u_vel - domainparticles[j].u_vel;
                dvy = domainparticles[i].v_vel - domainparticles[j].v_vel;
                dx = domainparticles[i].x_coordinate - domainparticles[j].x_coordinate;
                dy = domainparticles[i].y_coordinate - domainparticles[j].y_coordinate;
                sumx += ((domainparticles[j].mass)*(dvx*domainparticles[i].dwdx[j])*dx)/(domainparticles[j].rho*domainparticles[i].neighbour_distance[j]*domainparticles[i].neighbour_distance[j] + 0.01*h*h); //added 0.01h*h to avoid division by 0
                sumy += ((domainparticles[j].mass)*(dvy*domainparticles[i].dwdy[j])*dy)/(domainparticles[j].rho*domainparticles[i].neighbour_distance[j]*domainparticles[i].neighbour_distance[j] + 0.01*h*h);
            }
        }
      //  sumx = (2*0.000001*sumx)/domainparticles[i].rho; //viscosity added
      // sumy = (2*0.000001*sumy)/domainparticles[i].rho; 
    
         sumx = (2*0.0000008*sumx); //viscosity added. term 2 completed.
          sumy = (2*0.0000008*sumy); 
    
        domainparticles[i].dvxdt += sumx;
        domainparticles[i].dvydt += sumy; //gravity  removed 
       // sumx = 0;
       // sumy = 0;
    }
    for(int i = 0; i < num_inner_particles; i++){
        domainparticles[i].u_vel = domainparticles[i].u_vel + dt*(domainparticles[i].dvxdt);
        domainparticles[i].v_vel = domainparticles[i].v_vel + dt*(domainparticles[i].dvydt);
    }
}

void calculate_position(struct particle* domainparticles, int num_inner_particles, double dt)
{
    for(int i = 0; i < num_inner_particles; i++)
    {
            domainparticles[i].x_coordinate = domainparticles[i].x_coordinate + ( domainparticles[i].u_vel ) * dt  ;
            domainparticles[i].y_coordinate = domainparticles[i].y_coordinate + ( domainparticles[i].v_vel ) * dt  ;
        
    }

}




	



