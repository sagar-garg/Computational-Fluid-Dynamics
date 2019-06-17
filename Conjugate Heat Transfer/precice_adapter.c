#include <stdlib.h>
#include <stdio.h>
#include "precice/SolverInterfaceC.h"
#include "precice_adapter.h"
#include "helper.h"

//This function is  useful to save values of U,V,T which can be restored again in restore checkpoint function
//This type of checkpoints are required when we are solving implicit problem.
void write_checkpoint(double time, double **U, double **V, double **T, double *time_cp, double **U_cp, double **V_cp,
                      double **T_cp, int imax, int jmax)
{
    for (int i=0; i<=imax; i++)
    {
        for (int j=0; j<=jmax+1; j++)
        {
            U_cp[i][j]=U[i][j];
        }
    }

    for (int i=0; i<=imax+1; i++)
    {
        for (int j=0; j<=jmax; j++)
        {
            V_cp[i][j]=V[i][j];
        }
    }

    for (int i=0; i<=imax+1; i++)
    {
        for (int j=0; j<=jmax+1; j++)
        {
            T_cp[i][j]=T[i][j];
        }
    }
}


//This function is  useful to Restore values of U,V,T which has been stored by writecheckpoint function
//This type of checkpoints are required when we are solving implicit problem.
void restore_checkpoint(double *time, double **U, double **V, double **T, double time_cp, double **U_cp,double **V_cp,
                        double **T_cp, int imax, int jmax)
{
    for (int i=0; i<=imax; i++)
    {
        for (int j=0; j<=jmax+1; j++)
        {
            U[i][j]=U_cp[i][j];
        }
    }

    for (int i=0; i<=imax+1; i++)
    {
        for (int j=0; j<=jmax; j++)
        {
            V[i][j]=V_cp[i][j];
        }
    }

    for (int i=0; i<=imax+1; i++)
    {
        for (int j=0; j<=jmax+1; j++)
        {
            T[i][j]=T_cp[i][j];
        }
    }

}

void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature, int *vertexIDs,
                               int temperatureID, double **TEMP, int **FLAG,double Th,double Tc){
                            
                            int h = 0;
                            for(int i = 0; i <= imax+1; i++){
                                for(int j = 0; j <= jmax+1; j++){
                                    switch (FLAG[i][j])
                                    {
                                    case 544: // E 
                                        temperature[h] = TEMP[i+1][j];
                                        h++;
                                        break;
                                    case 288: // W 
                                        temperature[h] = TEMP[i-1][j];
                                        h++;
                                        break;
                                    case 160: // S 
                                        temperature[h] = TEMP[i][j-1];
                                        h++;
                                        break;
                                    case 96: // N
                                        temperature[h] = TEMP[i][j+1];
                                        h++;
                                        break;


                                    case 608: // EN 
                                        temperature[h] = (TEMP[i+1][j] + TEMP[i][j+1])/2;
                                        h++;
                                        break;
                                    case 670: // ES 
                                        temperature[h] = (TEMP[i+1][j] + TEMP[i][j-1])/2;
                                        h++;
                                        break;
                                    case 416: // WS
                                        temperature[h] = (TEMP[i-1][j] + TEMP[i][j-1])/2;
                                        h++;
                                        break;
                                    case 352: // WN
                                        temperature[h] = (TEMP[i-1][j] + TEMP[i][j+1])/2;
                                        h++;
                                        break;


                                    default:
                                        break;
                                    }      
                                    }
                            }
                            if(h==num_coupling_cells){
                                printf("num_coupling_cells is OK 3 \n");
                            }

                            //Convert dimentionless Temp. to Dimentioned Temp. before sending to PreCICE
                            //double Th = max_T(temperature, num_coupling_cells);
                            //double Th=10;
                            //double Tc = min_T(temperature, num_coupling_cells);
                            //double Tc=0;
                            double deltat = Th - Tc;
                            double Tinf = Tc;
                            for(int i = 0; i < num_coupling_cells; i++){
                                temperature[i] = (temperature[i]*deltat) + Tinf;

                            }
                            precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
                            
                            //free(temperature);

}
                
void set_coupling_boundary(int num_coupling_cells, int imax, int jmax, double dx, double dy, double *heatflux, double **TEMP, int **FLAG){                            
                            
                            int h = 0;
                            //int num_coupling_cells=1; //remove later
                            for(int i = 0; i <= imax+1; i++){
                                for(int j = 0; j <= jmax+1; j++){
                                    switch (FLAG[i][j])
                                    {
                                    
                                    case 544: // E 
                                        TEMP[i][j] = TEMP[i+1][j] + dx*heatflux[h];
                                        //TEMP[i][j] = TEMP[i+1][j] + dy*heatflux[h];
                                        h++;
                                        break;
                                    case 288: // W 
                                        TEMP[i][j] = TEMP[i-1][j] + dx*heatflux[h];
                                        //TEMP[i][j] = TEMP[i-1][j] + dy*heatflux[h];
                                        h++;
                                        break;
                                    case 160: // S 
                                        TEMP[i][j] = TEMP[i][j-1] + dy*heatflux[h];
                                        //TEMP[i][j] = TEMP[i][j-1] + dx*heatflux[h];
                                        h++;
                                        break;
                                    case 96: // N
                                        TEMP[i][j] = TEMP[i][j+1] + dy*heatflux[h];
                                        //TEMP[i][j] = TEMP[i][j+1] + dx*heatflux[h];
                                        h++;
                                        break;
                                    case 608: // EN 
                                        TEMP[i][j] = TEMP[i][j+1] + dy*heatflux[h]; //only considering north.
                                        //TEMP[i][j]=(TEMP[i+1][j] + dx*heatflux[h] + TEMP[i][j+1] + dy*heatflux[h])/2 ;
                                        h++;
                                        break;
                                    case 670: // ES 
                                        TEMP[i][j] = TEMP[i][j-1] + dy*heatflux[h];
                                        //TEMP[i][j]=(TEMP[i+1][j] + dx*heatflux[h] + TEMP[i][j-1] + dy*heatflux[h])/2;
                                        h++;
                                        break;
                                    case 416: // WS
                                        TEMP[i][j] = TEMP[i][j-1] + dx*heatflux[h];
                                        //TEMP[i-1][j] = (TEMP[i-1][j] + dx*heatflux[h] + TEMP[i][j-1] + dy*heatflux[h])/2;
                                        h++;
                                        break;
                                    case 352: // WN
                                        TEMP[i][j] = TEMP[i][j+1] + dx*heatflux[h];
                                        //TEMP[i][j] = (TEMP[i-1][j] + dx*heatflux[h] + TEMP[i][j+1] + dy*heatflux[h])/2;
                                        h++;
                                        break;

                                    default:
                                        break;
                                    }      
                                    }
                            }
                            if(h==num_coupling_cells){
                                printf("num_coupling_cells is OK 2 \n");
                            }
}




int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells, int meshID, int **FLAG)
{
    int dim= precicec_getDimensions();
                            //int dim=3;              
                            int h = 0;
                            
                            double v[num_coupling_cells*dim];
                            double *vertices=v;
                            int* vertexIDs=(int *) malloc(sizeof(int)*num_coupling_cells*dim); 
                            
                            for(int i = 0; i <= imax+1; i++){
                                for(int j = 0; j <= jmax+1; j++){
                                    switch (FLAG[i][j])
                                    {
                                    case 544: // E 
                                        vertices[h] = x_origin + i*dx;
//                                        printf("i= %d, j=%d \n", i, j);
                                        vertices[h+1]=y_origin + (j-0.5)*dy;
                                        vertices[h+2]=0;
                                        h=h+3;
                                        break;
                                    case 288: // W 
                                        vertices[h] = x_origin + i*dx;
  //                                      printf("i= %d, j=%d \n", i, j);
                                        vertices[h+1]=y_origin + (j-0.5)*dy;
                                        vertices[h+2]=0;
                                        h=h+3;
                                        break;
                                    case 160: // S 
                                        vertices[h] = x_origin + (i-0.5)*dx;
    //                                    printf("i= %d, j=%d \n", i, j);
                                        vertices[h+1]=y_origin + (j)*dy;
                                        vertices[h+2]=0;
                                        h=h+3;
                                        break;
                                    case 96: // N
                                        vertices[h] = x_origin + (i-0.5)*dx;
      //                                  printf("i= %d, j=%d \n", i, j);
                                        vertices[h+1]=y_origin + (j)*dy;
                                        vertices[h+2]=0;
                                        h=h+3;
                                        break;

                                    case 608: // EN   //only considering same as north
                                        vertices[h] = x_origin + (i-0.5)*dx;
      //                                  printf("i= %d, j=%d \n", i, j);
                                        vertices[h+1]=y_origin + (j)*dy;
                                        vertices[h+2]=0;
                                        h=h+3;
                                        break;
                                    case 670: // ES 
                                        vertices[h] = x_origin + (i-0.5)*dx;
    //                                    printf("i= %d, j=%d \n", i, j);
                                        vertices[h+1]=y_origin + (j)*dy;
                                        vertices[h+2]=0;
                                        h=h+3;
                                        break;
                                    case 416: // WS
                                        vertices[h] = x_origin + (i-0.5)*dx;
      //                                  printf("i= %d, j=%d \n", i, j);
                                        vertices[h+1]=y_origin + (j)*dy;
                                        vertices[h+2]=0;
                                        h=h+3;
                                        break;
                                    case 352: // WN
                                        vertices[h] = x_origin + (i-0.5)*dx;
      //                                  printf("i= %d, j=%d \n", i, j);
                                        vertices[h+1]=y_origin + (j)*dy;
                                        vertices[h+2]=0;
                                        h=h+3;
                                        break;


                                    default:
                                        break;
                                    }      
                                    }
                            }

                            if(h==(num_coupling_cells*dim)){
                                printf("num_coupling_cells is OK1 \n");
                            }
/*
                            for (int k=0 ; k<(num_coupling_cells*dim); k++)
                            {
                                printf("%lf ",vertices[k]);
                            }
                            printf("\n");
  */                          
                            precicec_setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
                            //free(vertices);                        
return(vertexIDs);
}

/*
int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells, int meshID, int **FLAG){
int dim = precice_getDimensions();
   // int dimension=3;
// int num_coupling_cells=30;
  //    int imax=15;
  //int jmax=15;

    //double v[90];
  
  //  double iscoupled[imax+2][jmax+2];
int** iscoupled=NULL;
  iscoupled = imatrix(0, (imax+2)-1, 0, (jmax+2)-1); //Initiation of Flag matrix.

  // double x_origin=1;
   // double y_origin=1;
  // double dx=1;
  // double dy=1;

  int i_sw, j_sw, i_se, j_se, i_ne, j_ne, i_nw, j_nw;

  for (int i=0; i<=imax+1; i++){
        for (int j=0; j<=jmax+1; j++){
            iscoupled[i][j]=0;
        }
        }

   for (int i=0; i<=imax+1; i++){
        for (int j=0; j<=jmax+1; j++){
            if (flag==)
            iscoupled[i][j]=0;
        }
        } 


/*
          for (int i=5; i<=11; i++){
              iscoupled[i][0]=1;
              iscoupled[i][7]=1;
          }

          for (int j=0; j<=7; j++){
              iscoupled[5][j]=1;
              iscoupled[11][j]=1;
          }

     for (int i=0; i<=imax+1; i++){
        for (int j=0; j<=jmax+1; j++){
            printf("%lf ", iscoupled[i][j]);

        }
        printf("\n");
            
        }     
*/

/*
int k=0;
//bottom
    for (int i=0; i<=imax; i++){
        if (iscoupled[i][0]==1){
            v[k]=k;
            v[k+1]=k+1;
            v[k+2]=k+2;
            k=k+3;
            i_se=i;
            j_se=0;
            printf("i_se =%d ", i_se);
             printf("j_se= %d ", j_se);
             printf("k= %d\n", k);     
        }
    }

//right
    for (int j=0; j<=jmax; j++){
        if (iscoupled[i_se][j]==1){
            v[k]=k;
            v[k+1]=k+1;
            v[k+2]=k+2;
            k=k+3;
            j_ne=j;
            i_ne=i_se;
            printf("i_ne =%d ", i_ne);
             printf("j_ne= %d ", j_ne);
                  printf("k= %d\n", k);
        
        }
    }
    

//top
    for (int i=i_ne; i>=0; i--){
        if (iscoupled[i][j_ne]==1){
            v[k]=k;
            v[k+1]=k+1;
            v[k+2]=k+2;
            k=k+3;
            i_nw=i;
            j_nw=j_ne;
            printf("i_nw =%d ", i_nw);
             printf("j_nw= %d ", j_nw);
                  printf("k= %d\n", k);
        
        }
    }

//left
    for (int j=j_nw; j>=0; j--){
        if (iscoupled[i_nw][j]==1){
            v[k]=k;
            v[k+1]=k+1;
            v[k+2]=k+2;
            k=k+3;
            i_sw=i_nw;
            j_sw=j_nw;
            printf("i_nw =%d ", i_sw);
             printf("j_nw= %d ", j_sw);
                  printf("k= %d\n", k);
        
        }
    }

}                           



}
*/                            
