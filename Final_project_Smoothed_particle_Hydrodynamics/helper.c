#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>
#include "helper.h"

/* ----------------------------------------------------------------------- */
/*                         local auxiliary functions                       */
/* ----------------------------------------------------------------------- */

clock_t last_timer_reset;

int min_int( const int n1, const int n2 )
{
    if( n1 < n2 ) return n1;
    return n2;
}

double fmin( double a, double b)
{
    if( a < b ) return a;
    return b;
}

double fmax( double a, double b)
{
    if( a > b ) return a;
    return b;
}

/* ----------------------------------------------------------------------- */


//Function to calculate the distance between two particles 
double calculate_distance(double x1, double y1, double x2, double y2)
{
  double dist = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
  return dist; 
}

//kernel function
double kernel (double dist, double h)
{
  double W=0;
 // double h=radius/k_scale;  //k_scale=2
  double factor=15.0/(7.0*(3.14)*h*h);
  double q=dist/h;
 // double l1=(2.0/3.0) - (q*q) + 0.5*q*q*q;
 // double l2=(1.0/6.0)*(2-q)*(2-q)*(2-q);
  //printf( "dist= %lf, factor= %lf, q= %lf, h=%lf, l1=%lf, l2-%lf \n", dist, factor, q, h, l1, l2);

  if (q<=1.0)
  {
    W=factor*( (2.0/3.0) - (q*q) + 0.5*q*q*q);
  }
  else if (q<=2.0)
  {  
    W=factor*(1.0/6.0)*(2.0-q)*(2.0-q)*(2.0-q);
  }
  else
  {
    W=0;
  }

return W;
}


//gradient of kernel function 
double kernel_grad (double dist, double h)
{
  double dwdq=0;
 // double h=radius/k_scale;  //k_scale=2
  double factor=15.0/(7.0*(3.14)*h*h);
  double q=dist/h;


  if (q<=1)
  {
    dwdq=factor*(-2.0*q+ 1.5*q*q)/(h);
  }
  else if (q<=2)
  {  
    dwdq=factor*(3.0*(2.0-q)*(2.0-q)*(-1.0/h));
  }
  else
  {
     dwdq=0;
  }
return dwdq;
}

double kernel_grad_x(double dist, double h, double x1, double x2)
{
  double dwdx=0;
 // double h=radius/k_scale;  //k_scale=2
  double factor=15.0/(7.0*(3.14)*h*h);
  double q = dist/h;
  double dx = x1-x2;


  if (q>0 && q<=1)
  {
  //dwdx = factor * (-2.0 + (1.5 * q))* dx;
  

  }
  else if (q>1 && q<=2)
  {  
    dwdx = -factor * 0.5* (2.0 - q) * (2.0 - q) * (dx/q);
  }
  else
  {
     dwdx=0;
  }
return dwdx;
}



void set_coordinates(struct particle* domainparticles, double x_len, double y_len, int num_boundary_particles, int num_inner_particles)
{
  int div = sqrt(num_inner_particles);  //40
  double hx = x_len/div; // spacing in x-direction 0.25
  double hy = y_len/div; // spacing in y-direction 0.25
  int count = 0;
  for(int j = 0; j < div; j++){
    for(int i = 0; i< div; i++){
      domainparticles[count].boundary = 0; // set non boundary particles
      domainparticles[count].particle_id = count; // set particle id e.g. = 0, 1, 2, 3 ...
      domainparticles[count].x_coordinate = i*hx + (hx/2.0); // for n = 3 (x,y) = (3.5*hx, hy)
      domainparticles[count].y_coordinate = j*hy + (hy/2.0);
      count = count + 1;
    }
  }


  int boundary_id = -1;
  for(int j = 0; j < div; j++){
      domainparticles[count].boundary = 1; // set boundary particles
      domainparticles[count].particle_id = boundary_id; // set particle id e.g. = -1, -2, -3 ...
      domainparticles[count].x_coordinate = j*hx + (hx/2.0); 
      domainparticles[count].y_coordinate = -(hy/2.0);
      boundary_id = boundary_id - 1;
      count = count + 1;
  }


  for(int j = 0; j < div; j++){
      domainparticles[count].boundary = 1; // set boundary particles
      domainparticles[count].particle_id = boundary_id; // set particle id e.g. = -1, -2, -3 ...
      domainparticles[count].x_coordinate = div*hx + (hx/2.0); 
      domainparticles[count].y_coordinate = j*hy + (hy/2.0);
      boundary_id = boundary_id - 1;
      count = count + 1;
  }


  for(int j = 0; j < div; j++){
      domainparticles[count].boundary = 1; // set boundary particles
      domainparticles[count].particle_id = boundary_id; // set particle id e.g. = -1, -2, -3 ...
      domainparticles[count].x_coordinate = (div-j-1)*hx + (hx/2.0); 
      domainparticles[count].y_coordinate = div*hy + (hy/2.0);
      boundary_id = boundary_id - 1;
      count = count + 1;
  }


  for(int j = 0; j < div; j++){
      domainparticles[count].boundary = 1; // set boundary particles
      domainparticles[count].particle_id = boundary_id; // set particle id e.g. = -1, -2, -3 ...
      domainparticles[count].x_coordinate = -(hx/2.0); 
      domainparticles[count].y_coordinate = (div-j-1)*hy + (hy/2.0);
      boundary_id = boundary_id - 1;
      count = count + 1;
  }
  
}


void reflect(struct particle* domainparticles, int num_inner_particles, double x_len, double y_len, double U_in)
{
  //int div = sqrt(num_inner_particles);  //40
  //double hx = x_len/div; // spacing in x-direction 0.25
  //double hy = y_len/div; // spacing in y-direction 0.25
  //int count=div;


  for (int i=0; i<num_inner_particles; i++)
  {
    double extra_dist=0;
// Coefficient of resitiution
//const double DAMP = 1;

//if (*u==0)
//return;
if (domainparticles[i].x_coordinate>x_len){
extra_dist=domainparticles[i].x_coordinate-x_len;
domainparticles[i].x_coordinate=x_len-extra_dist;
domainparticles[i].u_vel=-domainparticles[i].u_vel;
}

if (domainparticles[i].x_coordinate<0){
extra_dist=0-domainparticles[i].x_coordinate;
domainparticles[i].x_coordinate=0+extra_dist;
domainparticles[i].u_vel=-domainparticles[i].u_vel;
}

if (domainparticles[i].y_coordinate>y_len){
extra_dist=domainparticles[i].y_coordinate-y_len;
domainparticles[i].y_coordinate=y_len-extra_dist;
domainparticles[i].v_vel=-domainparticles[i].v_vel; //check
//domainparticles[i].u_vel=U_in;
}

if (domainparticles[i].y_coordinate<0){
extra_dist=0-domainparticles[i].y_coordinate;
domainparticles[i].y_coordinate=0+extra_dist;
domainparticles[i].v_vel=-domainparticles[i].v_vel;
}


  }
}


/* with damping
void reflect(struct particle* domainparticles, int num_inner_particles, double x_len, double y_len, double U_in)
{
  //int div = sqrt(num_inner_particles);  //40
  //double hx = x_len/div; // spacing in x-direction 0.25
  //double hy = y_len/div; // spacing in y-direction 0.25
  //int count=div;


  for (int i=0; i<num_inner_particles; i++)
  {
    double extra_dist=0;
// Coefficient of resitiution
//const double DAMP = 1;

//if (*u==0)
//return;
if (domainparticles[i].x_coordinate>x_len){
extra_dist=domainparticles[i].x_coordinate-x_len;
domainparticles[i].x_coordinate=x_len-0.75*extra_dist;
domainparticles[i].u_vel=-0.75*domainparticles[i].u_vel;
}

if (domainparticles[i].x_coordinate<0){
extra_dist=0-domainparticles[i].x_coordinate;
domainparticles[i].x_coordinate=0+0.75*extra_dist;
domainparticles[i].u_vel=-0.75*domainparticles[i].u_vel;
}

if (domainparticles[i].y_coordinate>y_len){
extra_dist=domainparticles[i].y_coordinate-y_len;
domainparticles[i].y_coordinate=y_len-0.75*extra_dist;
domainparticles[i].v_vel=-0.75*domainparticles[i].v_vel; //check
//domainparticles[i].u_vel=U_in;
}

if (domainparticles[i].y_coordinate<0){
extra_dist=0-domainparticles[i].y_coordinate;
domainparticles[i].y_coordinate=0+0.75*extra_dist;
domainparticles[i].v_vel=-0.75*domainparticles[i].v_vel;
}


  }
}
*/


//function to check the convergence
void check_convergence(struct particle* domainparticles, double* max_error, int num_inner_particles, int num_boundary_particle, double* old_temperature, double* old_u_velocity, double* old_v_velocity,int temp)
{
  double curr_error;
  
  if(temp==1)
  {
    *max_error = fabs(old_temperature[1]-(domainparticles[1].temperature));
    //printf("\nmax error = %lf",max_error);
    
    for (int i=2; i<num_inner_particles; i++)
    {
      curr_error = fabs(old_temperature[i]-(domainparticles[i].temperature));
      //printf("\nCurre error = %lf",curr_error);
      if( curr_error > *max_error )
      {
        *max_error = curr_error ;
      }      
    }
  }

  else 
  {
    *max_error = fmax(    fabs(old_u_velocity[1]-(domainparticles[1].u_vel)) ,   fabs(old_v_velocity[1]-(domainparticles[1].v_vel))   );
    //printf("\nmax error = %lf",max_error);
    
    for (int i=2; i<num_inner_particles; i++)
    {
      curr_error = fmax(    fabs(old_u_velocity[i]-(domainparticles[i].u_vel)) ,   fabs(old_v_velocity[i]-(domainparticles[i].v_vel))   );
      //printf("\nCurre error = %lf",curr_error);
      if( curr_error > *max_error )
      {
        *max_error = curr_error ;
      }      
    }
  }


}

//Cleaning all the dynamic memory 
void cleaning(struct particle* domainparticles, int num_particles, double* old_temperature, double* old_v_velocity, double* old_u_velocity, double* temp_laplacian)
{
for (int i=0; i<=num_particles-1; i++)
  {
    free(domainparticles[i].neighbour_distance);
    free(domainparticles[i].weights); 
    free(domainparticles[i].dwdx);
    free(domainparticles[i].dwdy);
  }

  free(domainparticles);
  free(old_temperature);
  free(old_u_velocity);
  free(old_v_velocity);
  free(temp_laplacian);
}



/* ----------------------------------------------------------------------- */
/*                             read datafile                               */
/* ----------------------------------------------------------------------- */

void errhandler( int nLine, const char *szFile, const char *szString )
{
    int err = errno;

    fprintf( ERROUT, "%s:%d Error : %s", szFile, nLine, szString );
    fprintf( ERROUT, "\n" );
    
    /* if an error within the c-library occured, an error code can be   */
    /* found in the global variable err                                 */
    if( err != 0 )
    {
	fprintf( ERROUT, "C-Lib   errno    = %d\n", err);
	fprintf( ERROUT, "C-Lib   strerror = %s\n", strerror( err ) );
    }
    exit(1);
}


/*  for comfort */
#define READ_ERROR(szMessage, szVarName, szFileName, nLine) \
  { char szTmp[80]; \
    if( nLine ) \
	sprintf( szTmp, " %s  File: %s   Variable: %s  Line: %d", szMessage, szFileName, szVarName, nLine ); \
    else \
	sprintf( szTmp, " %s  File: %s   Variable: %s ", szMessage, szFileName, szVarName); \
    ERROR( szTmp ); \
  }
    



/* --------------------------------------------------------------------------*/
/* The function searches the datafile fh for the line defining the variable  */
/* szVarName and returns the respctive string including the value of the     */
/* variable. If there's no appropriate line within the datafile, the program */
/* stops with an error messsage.                                             */
/* ATTENTION: The pointer returned refers to a static variable within the    */
/* function. To maintain the string over several program calls, it has to be */
/* copied!!!                                                                 */
/*                                                                           */
char* find_string( const char* szFileName, const char *szVarName )
{ 
    int nLine = 0;
    int i;
    FILE *fh = NULL;
    
    static char szBuffer[MAX_LINE_LENGTH];	/* containes the line read  */
                                               /* from the datafile        */

    char* szLine = szBuffer;
    char* szValue = NULL;
    char* szName = NULL;

    /* open file */
    fh = fopen( szFileName, "rt" );
    if( fh == 0 ) 
	READ_ERROR("Could not open file", szVarName, szFileName, 0);

    /* searching */
    while( ! feof(fh) )
    {
	fgets( szLine, MAX_LINE_LENGTH, fh );
	++nLine;

	/* remove comments */
	for( i = 0; i < strlen(szLine); i++)
	    if( szLine[i] == '#' )
	    {
		szLine[i] = '\0'; /* Stringende setzen */
		break;
	    }

	/* remove empty lines */
	while( isspace( (int)*szLine ) && *szLine) ++szLine;
	if( strlen( szLine ) == 0) continue; 

	/* now, the name can be extracted */
	szName = szLine;
	szValue = szLine;
	while( (isalnum( (int)*szValue ) || *szValue == '_') && *szValue) ++szValue;
	
	/* is the value for the respective name missing? */
	if( *szValue == '\n' || strlen( szValue) == 0)  
	    READ_ERROR("wrong format", szName, szFileName, nLine);
	
	*szValue = 0;		/* complete szName! at the right place */
	++szValue;
        
	/* read next line if the correct name wasn't found */
	if( strcmp( szVarName, szName)) continue;

	/* remove all leading blnkets and tabs from the value string  */
	while( isspace( (int)*szValue) ) ++szValue;
	if( *szValue == '\n' || strlen( szValue) == 0)  
	    READ_ERROR("wrong format", szName, szFileName, nLine);
	
	fclose(fh);
	return szValue;
    }  
   
    READ_ERROR("variable not found", szVarName, szFileName, nLine);
    
    return NULL;		/* dummy to satisfy the compiler  */
} 

void read_string( const char* szFileName, const char* szVarName, char*   pVariable)
{
    char* szValue = NULL;	/* string containg the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as variable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%s", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName,0);

    printf( "File: %s\t\t%s%s= %s\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              pVariable );
}

void read_int( const char* szFileName, const char* szVarName, int* pVariable)
{
    char* szValue = NULL;	/* string containing the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%d", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName, 0);

    printf( "File: %s\t\t%s%s= %d\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              *pVariable );
}

void read_double( const char* szFileName, const char* szVarName, double* pVariable)
{
    char* szValue = NULL;	/* String mit dem eingelesenen Variablenwert */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%lf", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName, 0);

    printf( "File: %s\t\t%s%s= %f\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              *pVariable );
}
