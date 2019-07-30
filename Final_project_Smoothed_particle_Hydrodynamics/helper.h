#ifndef __HELPER_H__
#define __HELPER_H__

/* includefiles */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

//Defining struct Particle 
struct particle 
{
  int particle_id;
  int boundary; //Flag=1 if this is boundary point
  double x_coordinate; 
  double y_coordinate;
  double u_vel;
  double v_vel;
  double pressure;
  double rho;
  double mass;
  double temperature;
  double *neighbour_distance; //array conatains the distance of each neighbouring particles
  double *weights; //array contains the weights of each neighbouring particles after applying Kernal function 
  double *dwdq;
  double *dwdx;
  double *dwdy;
  //double vx;
  //double vy;
  double drhodt;
  double dvxdt;
  double dvydt;
};

//Function to calculate the distance between two particles 
double calculate_distance(double x1, double x2, double y1, double y2);

//kernel function
double kernel (double q, double h);

//gradient of kernel function 
double kernel_grad (double q, double h);

double kernel_grad_x (double dist, double h, double x1, double x2);

//function to initialise coordinats of all particles
void set_coordinates(struct particle* domainparticles, double x_len, double y_len, int num_boundary_particles, int num_non_boundary_particles);

//function to set boundary for 
void reflect(struct particle* domainparticles, int num_inner_particles, double x_len, double y_len, double U_in);

//function to check the convergence
void check_convergence(struct particle* domainparticles, double* max_error, int num_inner_particles, int num_boundary_particle, double* old_temperature, double* old_u_velocity, double* old_v_velocity,int temp);

//Cleaning all the dynamic memory 
void cleaning(struct particle* domainparticles, int num_particles, double* old_temperature, double* old_v_velocity, double* old_u_velocity, double* temp_laplacian);

#define MAX_LINE_LENGTH 1024

/**
 * Stores the last timer value 
 */
extern clock_t last_timer_reset; 

/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
#define ERROR(s)    errhandler( __LINE__, __FILE__, s)

/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
#define ERROUT stdout

/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
void  errhandler( int nLine, const char *szFile, const char *szString );


/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_INT( szFileName, VarName)    read_int   ( szFileName, #VarName, &(VarName) ) 

/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_DOUBLE( szFileName, VarName) read_double( szFileName, #VarName, &(VarName) )

/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_STRING( szFileName, VarName) read_string( szFileName, #VarName,  (VarName) )

void read_string( const char* szFilename, const char* szName, char*  sValue);
void read_int   ( const char* szFilename, const char* szName, int*    nValue);
void read_double( const char* szFilename, const char* szName, double*  Value);


#endif

