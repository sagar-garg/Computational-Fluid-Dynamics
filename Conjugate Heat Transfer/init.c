#include "helper.h"
#include "init.h"

int read_parameters(const char *szFileName,    /* name of the file */
                    char *geometry,            /* name of the pgm geometry file */
                    char *precice_config,      /* name of the precice Configuration file */
                    char *participant_name,    /* name of the participant typically fluid*/
                    char *mesh_name,           /* Mesh name typically fluid_mesh */
                    char *read_data_name,      /* Parameter which is to be read from PreCICE typcally heat flux */
                    char *write_data_name,     /* Parameter which is to be pass to PreCICE tpically Temperature */
                    double *Re,                /* reynolds number   */
                    double *Pr,                /* Prandle number */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *U_in,              /* U_velocity inflow */
                    double *V_in,              /* V_velocity inflow */
                    double *PI,                /* pressure */
                    double *TI,                /* Temperature */
                    double *Te,                /* East side temperature */
                    double *Tw,                /* West side temperature */
                    double *Tn,                /* North side temperature */
                    double *Ts,                /* South side temperature */   
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *x_origin,          /* Xorigin value for the Coupling domain */
                    double *y_origin,          /* Yorigin value for the Coupling domain */
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    double *beta,              /* Coefficient of thermal expansion */
                    int  *itermax,             /* max. number of iterations- for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
		                double *dt_value)          /* time for output */
{
   READ_STRING( szFileName, geometry );
   READ_STRING( szFileName, precice_config );
   READ_STRING( szFileName, participant_name );
   READ_STRING( szFileName, mesh_name );
   READ_STRING( szFileName, read_data_name );
   READ_STRING( szFileName, write_data_name );

   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );
   READ_DOUBLE( szFileName, *x_origin );
   READ_DOUBLE( szFileName, *y_origin );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *Pr    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );   

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );
   READ_DOUBLE( szFileName, *beta    );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *U_in );
   READ_DOUBLE( szFileName, *V_in );

   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *TI );
   READ_DOUBLE( szFileName, *Te );
   READ_DOUBLE( szFileName, *Tw );
   READ_DOUBLE( szFileName, *Tn );
   READ_DOUBLE( szFileName, *Ts );

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;
}

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
  double **T_cp)
{
  init_matrix( P , 0 , imax+1 , 0 , jmax+1, PI ); /*Initiation of Pressure value matrix with PI value*/
    
  init_matrix( U, 0 , imax , 0 , jmax+1, UI);  /*Initiation of X_Velocity value matrix with UI value*/  

  init_matrix( V, 0 , imax+1 , 0 , jmax, VI);  /*Initiation of Y_Velocity value matrix with VI value*/
  
  init_matrix( T, 0, imax+1, 0 , jmax+1, TI );  /*Initiation of Temperature value matrix with TI value*/

init_matrix( U_cp, 0 , imax , 0 , jmax+1, 0);  /*Initiation of X_Velocity value matrix with UI value*/  

  init_matrix( V_cp, 0 , imax+1 , 0 , jmax, 0);  /*Initiation of Y_Velocity value matrix with VI value*/
  
  init_matrix( T_cp, 0, imax+1, 0 , jmax+1, 0 );  /*Initiation of Temperature value matrix with TI value*/

}

//This is file init_flag() which is used for making the int**flag matrix 
//after reading given Pic (from pgm file)

void init_flag(const char *szFileName,char* geometry, int imax, int jmax,int** flag,int* count_num_fluid,int* count_coupling_num)
{
  int** pic=NULL;
  int xsize=imax+2;
  int ysize=jmax+2;
  
  int** dummy=NULL; //This additional matrix is made to set the next set of directional bands
  dummy = imatrix(0, xsize-1, 0, ysize-1);
  pic=read_pgm((const char*)geometry); //Reading geometry.pgm file

  int i=0;
  while(i<xsize)
  {
    int j=0;
    while(j<ysize)
    {
      switch (pic[i][j])
      {
      case 0://0 represents No-slip in the pgm file
        flag[i][j]=1<<1;
        dummy[i][j]=0;//No fluid it will remain=0
        break;
      case 1://1 represents Free-slip in the pgm file 
        flag[i][j]=1<<2;
        dummy[i][j]=0;//No fluid it will remain=0
        break;
      case 2://2 represents Noflow
        flag[i][j]=1<<3;
        dummy[i][j]=0;//No fluid it will remain=0
        break;
      case 3://3 represents Inflow 
        flag[i][j]=1<<4;
        dummy[i][j]=0;//No fluid it will remain=0
        break;
      case 4://4 represents coupling cell
        flag[i][j]=1<<5;
        dummy[i][j]=0;//No fluid it will remain=0
        (*count_coupling_num)++;
        break;
      default://6 represents Fluid in the pgm file 
        flag[i][j]=1;//000000001 As it is Fluid
        dummy[i][j]=1;//It is for fluid.
        (*count_num_fluid)++;
        break;
      }
      j++;
    }
    i++;
  }

// Code for setting directional Flags North, South , East and West
//(Old)New = Binary representation - Meaning in the Geometry file
//(32)64 ``=     100 0000 - Fluid in North direction
//(64)128 `=    1000 0000 - Fluid in South direction
//(128)256 =  1 0000 0000 - Fluid in West direction
//(256)512 = 10 0000 0000 - Fluid in East direction


//Setting values of all Inner boundary cells 
for(int i=1;i<xsize-1;i++)
{
  for(int j=1;j<ysize-1;j++)
  {
    if(flag[i][j]!=1)//Only changes on Non-Fluid cells ie. on obstacle cells
    {
      flag[i][j]=  (flag[i][j])  |  ((dummy[i+1][j])*(512))  |  ((dummy[i-1][j])*(256))  |  
                 ((dummy[i][j+1])*(64))  |  ((dummy[i][j-1])*(128));
    }
  }
}


//Setting values of all Left and right boundary cells 
for(int j=1;j<ysize-1;j++)
{
  if(flag[0][j]!=1)//Only changes on Non-Fluid cells ie. on Left boundary obstacle cells
    {
      flag[0][j]=  (flag[0][j])  |  ((dummy[0+1][j])*(512))  |     
                 ((dummy[0][j+1])*(64))  |  ((dummy[0][j-1])*(128));
    }
  if(flag[xsize-1][j]!=1)//Only changes on Non-Fluid cells ie. on right boundary obstacle cells
    {
      flag[xsize-1][j]=  (flag[xsize-1][j])   |  ((dummy[xsize-1-1][j])*(256))  |  
                 ((dummy[xsize-1][j+1])*(64))  |  ((dummy[xsize-1][j-1])*(128));
    }
}


//Setting values of all Bottom and TOP boundary cells 
for(int i=1;i<xsize-1;i++)
{
  if(flag[i][0]!=1)//Only changes on Non-Fluid cells ie. on Bottom boundary obstacle cells
    {
      flag[i][0]=  (flag[i][0])  |  ((dummy[i+1][0])*(512))  |  ((dummy[i-1][0])*(256))  |  
                 ((dummy[i][0+1])*(64));
      
    }
  if(flag[i][ysize-1]!=1)//Only changes on Non-Fluid cells ie. on TOP boundary obstacle cells
    {
      flag[i][ysize-1]=  (flag[i][ysize-1])  |  ((dummy[i+1][ysize-1])*(512))  |  ((dummy[i-1][ysize-1])*(256))   
                 |  ((dummy[i][ysize-1-1])*(128));  
    }
}

//Setting values at corner cell
flag[0][0]=(flag[0][0])  |  ((dummy[1][0])*(512))  |  ((dummy[0][1])*(64));//Bottom left
flag[0][ysize-1]=(flag[0][ysize-1])  |  ((dummy[1][ysize-1])*(512))   |  ((dummy[0][ysize-1-1])*(128));//Top Left
flag[xsize-1][0]=(flag[xsize-1][0])  |  ((dummy[xsize-1-1][0])*(256))  |  ((dummy[xsize-1][0+1])*(64));//Bottom Right
flag[xsize-1][ysize-1]=(flag[xsize-1][ysize-1])  |  ((dummy[xsize-1-1][ysize-1])*(256))  |  ((dummy[xsize-1][ysize-1-1])*(128));//Top Right

}
