#ifndef __INIT_H_
#define __INIT_H_

/**
 * This operation initializes all the local variables reading a configuration
 * file. For every variable a macro like READ_INT() is called passing it the
 * input filename and the variable to be written to. This macro calls
 * an operation read_int() augmenting the parameter set with the name of the
 * variable to be read. The read_int() operation parses the input file, extracts
 * the value of the variable, sets the variable and finally prints some debug
 * information. This is possible as the macro adds the name of the variable to
 * be set. All the helper operations can be found within helper.h and helper.c.
 *
 * @param szFileName char pointer to the filename
 * @param geometry   name of the pgm geometry file 
 * @param precice_config      name of the precice Configuration file 
 * @param participant_name     name of the participant typically fluid
 * @param mesh_name,           Mesh name typically fluid_mesh 
 * @param read_data_name,      Parameter which is to be read from PreCICE typcally heat flux 
 * @param write_data_name,     Parameter which is to be pass to PreCICE tpically Temperature 
 * @param Re         Reynolds number
 * @param UI         initial velocity in  x-direction - used by init_uvp()
 * @param VI         initial velocity y-direction - used by init_upv()
 * @param PI         initial pressure - used by init_upv()
 * @param TI         initial Temperature - used by init_upv()
 * @param Th         Hot side Boundary temperature
 * @param Tc         Cold side Boundary temperature
 * @param GX         gravitation x-direction
 * @param GY         gravitation y-direction
 * @param t_end      end time (not discrete in time steps)
 * @param xlength    domain length x-direction
 * @param ylength    domain lenght y-direction
 * @param x_origin   Xorigin value for the Coupling domain 
 * @param y_origin   Yorigin value for the Coupling domain
 * @param dt         time step length: dividing t_end by dt gives the number of
 *                   time steps to perform. Actually dt is determined by a
 *                   function, so manipulating this value within the 
 *                   configuration file should not affect the solution process
 *                   at all
 * @param dx         cell length x-direction
 * @param dy         cell length y-direction
 * @param imax       number of cells in x-direction
 * @param jmax       number of cells in Y-direction
 * @param alpha      uppwind-differencing-factor alpha
 * @param omg        relaxation factor omega
 * @param tau        safety parameter for time step calculation
 * @param beta       Thermal Expanson Coefficient
 * @param itermax    max. number of pressure iterations
 * @param eps        tolerance limit for pressure calculation
 * @param dt_value   time steps for output (after how many time steps one should
 *                   write into the output file)
 */
int read_parameters( 
  const char *szFileName,
  char *precice_config,      
  char *participant_name,    
  char *mesh_name,           
  char *read_data_name,      
  char *write_data_name,    
  char *geometry,
  double *Re,
  double *Pr,
  double *UI,
  double *VI,
  double *U_in,
  double *V_in,
  double *PI,
  double *TI,
  double *Te,
  double *Tw,
  double *Tn,
  double *Ts,
  double *GX,
  double *GY,
  double *t_end,
  double *xlength,
  double *ylength,
  double *x_origin,          
  double *y_origin,          
  double *dt,
  double *dx,
  double *dy,
  int  *imax,
  int  *jmax,
  double *alpha,
  double *omg,
  double *tau,
  double *beta,
  int  *itermax,
  double *eps,
  double *dt_value
);

/**
 * The arrays U,V,P and T are initialized to the constant values UI, VI , PI, TI on
 * the whole domain.
 */
void init_uvp(
  double UI,
  double VI,
  double PI,
  double TI,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P,
  double **T,
  double **U_cp,
    double **V_cp,
  double **T_cp

);


/**
 * This function set values to Flag matrix according to the boundary given 
 * in the pic matrix(taken from geometry.pgm)
 */

void init_flag(
  const char *szFileName,
  char* geometry, 
  int xsize, 
  int ysize,
  int** flag,
  int* count_num_fluid,
  int* count_coupling_num);

#endif

