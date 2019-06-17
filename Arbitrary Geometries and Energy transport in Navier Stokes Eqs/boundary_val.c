//  boundary_val.c

#include "boundary_val.h"
#include <stdio.h>
#include <assert.h>

//T_h is inflow, T_c is also inflow (Dirichlet). there is no outflow
//we assume all the wall is inflow (dirichlet), or insulated (adiabatic)
//in all questions related to Temp, we assume initial values of U, V are zero. 
//but all walls are no-slip or free-slip?

//spec boudnary val does not need Ui, VI

// in spec boundary val should we put i<imax or i<=imax

//for non-obstacle boundaries
void boundaryvalues(
                    int imax,
                    int jmax,
                    double **U,
                    double **V, int **flag, double UI, double VI, double U_in, double V_in)
{   
    //horizontal walls
    int i=1;
    while (i<=imax){
        
        //bottom boundary
        
        switch(flag[i][0])
        {
    case 258: case 130: case 66: case  34: case  290: case  322: case  194: case  162: //no slip
        //bottom boundary
        V[i][0]=0;
        U[i][0]=-U[i][1];
        break;
        
    case 260: case  132: case  68: case  36: case  292: case  324: case  196: case  164://free slip
        //bottom boundary
        V[i][0]=0;
        U[i][0]=U[i][1];
        break;
        
    case 264: case  136: case  72: case  40: case  296: case  328: case  200: case  168: //outflow. (should inflow is considered in spec_boundary_val)
        //bottom boundary
        U[i][0]=U[i][1];
        V[i][0]=V[i][1];
        break;
        
    case 272: case  144: case  80: case  48: case  304: case  336: case  208: case  176: //inflow
        //bottom
        U[i][0]=2*(U_in)- U[i][1];
        V[i][0]=V_in;
        break; //may need changes
    }
    i++;
}

//top boundary
    i=1;
    while (i<=imax){
        
    switch(flag[i][jmax+1])
    {
        case 258: case 130: case 66: case  34: case  290: case  322: case  194: case  162: //no slip
            //top boundary
            V[i][jmax]=0;
            U[i][jmax+1]=-U[i][jmax];
            break;
            
        case 260: case  132: case  68: case  36: case  292: case  324: case  196: case  164://free slip
            //top boundary
            V[i][jmax]=0;
            U[i][jmax+1]=U[i][jmax];
            break;
            
        case 264: case  136: case  72: case  40: case  296: case  328: case  200: case  168: //outflow. (should inflow is considered in spec_boundary_val)
            //top boundary
            U[i][jmax+1]=U[i][jmax];
            V[i][jmax]=V[i][jmax-1];
            break;
            
        case 272: case  144: case  80: case  48: case  304: case  336: case  208: case  176: //inflow
            //top
            U[i][jmax+1]=2*(U_in)-U[i][jmax];
            V[i][jmax]=(V_in);
            break; //may need changes
    }
    i++;
}



//for left boundary walls
int j=1;
while (j<=jmax) {
    switch(flag[0][j])
    {
        case 258: case 130: case 66: case  34: case  290: case  322: case  194: case  162: //no slip
            //left boundary
            U[0][j]=0;
            V[0][j]=-V[1][j];
            break;
            
        case 260: case  132: case  68: case  36: case  292: case  324: case  196: case  164://free slip
            //left boundary
            U[0][j]=0;
            V[0][j]=V[1][j];
            break;
            
        case 264: case  136: case  72: case  40: case  296: case  328: case  200: case  168: //outflow. (should inflow is considered in spec_boundary_val)
            //left boundary
            U[0][j]=U[1][j];
            V[0][j]=V[1][j];
            break;
            
        case 272: case  144: case  80: case  48: case  304: case  336: case  208: case  176: //inflow
            //left
            U[0][j]=U_in;
            V[0][j]=2*(V_in)-V[1][j];
            break; //may need changes
    }
    j++;
}

//for right boundary walls
j=1;
while (j<=jmax){
    switch(flag[imax+1][j])
    {
        case 258: case 130: case 66: case  34: case  290: case  322: case  194: case  162: //no slip
            //right boundary
            U[imax][j]=0;
            V[imax+1][j]=-V[imax][j];
            
            break;
            
        case 260: case  132: case  68: case  36: case  292: case  324: case  196: case  164://free slip
            //right boundary
            U[imax][j]=0;
            V[imax+1][j]=V[imax][j];
            break;
            
        case 264: case  136: case  72: case  40: case  296: case  328: case  200: case  168: //outflow.   //right boundary
            U[imax][j]=U[imax-1][j];
            V[imax+1][j]=V[imax][j];
            break;
            
        case 272: case  144: case  80: case  48: case  304: case  336: case  208: case  176: //inflow
            //right
            //U[imax-1][j]=U_in;
            U[imax][j]=U_in;
            //V[imax][j]=2*(V_in)-V[imax-1][j];
            V[imax+1][j]=2*(V_in)-V[imax][j];
            break; //may need changes
    }
    j++;
}

}

// add assert function
//should we add pressure here?
//for non-obstacle boundaries
//This is with Temperature
void boundaryvalues_t(
                    int imax,
                    int jmax,
                    double **U,
                    double **V,double **T, int **flag, double UI, double VI,double U_in, double V_in,double Te,double Tw, double Tn, double Ts)
{ // outflow, no-slip, free-slip.
    //inflow to be read in spec_boundary_val
    
    //horizontal walls
    int i=1;
    while (i<=imax){
        
        //bottom boundary
        if (Ts!=1000)
        {
            T[i][0]=2*Ts-T[i][1];//Dirichlet Temperature Boundary condition
        }
        else //Wall is adiabatic
        {
            T[i][0]=T[i][1]; //Neumann Temperature Boundary condition 
        }


        switch(flag[i][0])
        {
    case 258: case 130: case 66: case  34: case  290: case  322: case  194: case  162: //no slip
        //bottom boundary
        V[i][0]=0;
        U[i][0]=-U[i][1];
        break;
        
    case 260: case  132: case  68: case  36: case  292: case  324: case  196: case  164://free slip
        //bottom boundary
        V[i][0]=0;
        U[i][0]=U[i][1];
        break;
        
    case 264: case  136: case  72: case  40: case  296: case  328: case  200: case  168: //outflow. (should inflow is considered in spec_boundary_val)
        //bottom boundary
        U[i][0]=U[i][1];
        V[i][0]=V[i][1];
        break;
        
    case 272: case  144: case  80: case  48: case  304: case  336: case  208: case  176: //inflow
        //bottom
        U[i][0]=2*(U_in)- U[i][1];
        V[i][0]=(V_in);
        break; //may need changes
    }
    i++;
}

//top boundary
    i=1;
    while (i<=imax){
    
    if (Tn!=1000)
        {
            T[i][jmax+1]=2*Tn-T[i][jmax]; //Dirichlet Temperature Boundary condition
        }
        else //Wall is adiabatic
        {
            T[i][jmax+1]=T[i][jmax]; //Neumann Temperature Boundary condition 
        }

    switch(flag[i][jmax+1])
    {
        case 258: case 130: case 66: case  34: case  290: case  322: case  194: case  162: //no slip
            //top boundary
            V[i][jmax]=0;
            U[i][jmax+1]=-U[i][jmax];
            break;
            
        case 260: case  132: case  68: case  36: case  292: case  324: case  196: case  164://free slip
            //top boundary
            V[i][jmax]=0;
            U[i][jmax+1]=U[i][jmax];
            break;
            
        case 264: case  136: case  72: case  40: case  296: case  328: case  200: case  168: //outflow. (should inflow is considered in spec_boundary_val)
            //top boundary
            U[i][jmax+1]=U[i][jmax];
            V[i][jmax]=V[i][jmax-1];
            break;
            
        case 272: case  144: case  80: case  48: case  304: case  336: case  208: case  176: //inflow
            //top
            U[i][jmax+1]=2*(U_in)-U[i][jmax];
            V[i][jmax]=(V_in);
            break; //may need changes
    }
    i++;
}



//for left boundary walls
int j=1;
while (j<=jmax) {
    
    if (Tw!=1000)
        {
            T[0][j]=2*Tw-T[1][j]; //Dirichlet Temperature Boundary condition
        }
        else //Wall is adiabatic
        {
            T[0][j]=T[1][j]; //Neumann Temperature Boundary condition 
        }

    switch(flag[0][j])
    {
        case 258: case 130: case 66: case  34: case  290: case  322: case  194: case  162: //no slip
            //left boundary
            U[0][j]=0;
            V[0][j]=-V[1][j];
            break;
            
        case 260: case  132: case  68: case  36: case  292: case  324: case  196: case  164://free slip
            //left boundary
            V[0][j]=V[1][j];
            U[0][j]=0;
            break;
            
        case 264: case  136: case  72: case  40: case  296: case  328: case  200: case  168: //outflow. (should inflow is considered in spec_boundary_val)
            //left boundary
            U[0][j]=U[1][j];
            V[0][j]=V[1][j];
            break;
            
        case 272: case  144: case  80: case  48: case  304: case  336: case  208: case  176: //inflow
            //left
            U[0][j]=(U_in);
            V[0][j]=2*(V_in)-V[1][j];
            break; //may need changes
    }
    j++;
}

//for right boundary walls
j=1;
while (j<=jmax){
    
    if (Te!=1000)
        {
            T[imax+1][j]=2*Te-T[imax][j]; //Dirichlet Temperature Boundary condition
        }
        else //Wall is adiabatic
        {
            T[imax+1][j]=T[imax][j]; //Neumann Temperature Boundary condition 
        }

    
    
    switch(flag[imax+1][j])
    {
        case 258: case 130: case 66: case  34: case  290: case  322: case  194: case  162: //no slip
            //right boundary
            U[imax][j]=0;
            V[imax+1][j]=-V[imax][j];
            break;
            
        case 260: case  132: case  68: case  36: case  292: case  324: case  196: case  164://free slip
            //right boundary
            U[imax][j]=0;
            V[imax+1][j]=V[imax][j];
            break;
            
        case 264: case  136: case  72: case  40: case  296: case  328: case  200: case  168: //outflow.   //right boundary
            U[imax][j]=U[imax-1][j];
            V[imax+1][j]=V[imax][j];
            break;
            
        case 272: case  144: case  80: case  48: case  304: case  336: case  208: case  176: //inflow
            //right
            //U[imax-1][j]=U_in;
            U[imax][j]=U_in;
            //V[imax][j]=2*(V_in)-V[imax-1][j];
            V[imax+1][j]=2*(V_in)-V[imax][j];
            break; //may need changes
    }
    j++;
}

}

//for obstacle boundaries
void spec_boundary_val(int imax, int jmax, double **U, double **V, int **flag, double UI, double VI) {
    
    //TBD: include assert() function to give error if Boundary with 3 or 4 sharing walls present
    // assert (condition && "sorry dude! can't work with obstacles having opposite or 3 or 4 boundaries with fluid");
    
    /// run for all edge and corner cells only
    
    int i=1; //since only running for obstacles, can start form 1
    while (i<=imax)
    {
        int j=1;
        while (j<=jmax)
        {
            
            //if an edge . cant have edge and corner in the same switch statements as one may over-rule the other
            //if (flag==edge)
            //{
            switch (flag[i][j])
            {
                case 34: //BN
                    V[i][j]=0;
                    U[i-1][j]=-U[i-1][j+1];
                    U[i][j]=-U[i][j+1];
                    break;
                case 130: //bw
                    U[i-1][j]=0;
                    V[i][j-1]=-V[i-1][j-1];
                    V[i][j]=-V[i-1][j];
                    break;
                case 66: //bs
                    V[i][j-1]=0;
                    U[i-1][j]=-U[i-1][j-1];
                    U[i][j]=-U[i][j-1];
                    break;
                case 258: //bo
                    U[i][j]=0;
                    V[i][j]=-V[i+1][j];
                    V[i][j-1]=-V[i+1][j-1];
                    break;
                    //   default: //need to change this to show error
                    //       printf ("warning: not an edge. check again");
                    //       break;
                    //    }
                    //BN
                    //}
                    //else if (flag==corner)  //could also use else
                    //{
                    //  switch (flag[i][j])
                    //     {
                    
                    //corner cells
                case 290: //BN, BO
                    U[i][j]=0;
                    V[i][j]=0;
                    U[i-1][j]=-U[i-1][j+1];
                    V[i][j-1]=-V[i+1][j-1];
                    break;
                    
                case 162: //bn, bw
                    U[i-1][j]=0;
                    V[i][j]=0;
                    U[i][j]=-U[i][j+1];
                    V[i][j-1]=-V[i-1][j-1];
                    break;
                    
                case 194: //bw, bs
                    U[i-1][j]=0;
                    V[i][j-1]=0;
                    U[i][j]=-U[i][j-1];
                    V[i][j]=-V[i-1][j];
                    break;
                    
                case 322: //bs, bo
                    U[i][j]=0;
                    V[i][j-1]=0;
                    U[i-1][j]=-U[i-1][j-1];
                    V[i][j]=-V[i+1][j];
                    break;
                    
                    //      default: //need to change this to show error
                    //       printf ("warning: not a corner. check again");
                    //        break; //this is for fluids
            }
            
            //}
            //   else
            //  {
            //        break;// do nothing if fluid cell
            //   }
            j++;
        }
        i++;
    }
}


//With Temperature
//for obstacle boundaries
void spec_boundary_val_t(int imax, int jmax, double **U, double **V, int **flag, double **T) {
    //TBD: include assert() function to give error if Boundary with 3 or 4 sharing walls present
    // assert (condition && "sorry dude! can't work with obstacles having opposite or 3 or 4 boundaries with fluid");
    
    /// run for all edge and corner cells only
  
    int i=1; //since only running for obstacles, can start form 1
    while (i<=imax)// (< or <=   ???)
    {
        int j=1;
        while (j<=jmax)
        {
            switch (flag[i][j])
            {
                //edge cells
                case 34: //BN
                    V[i][j]=0;
                    U[i-1][j]=-U[i-1][j+1];
                    U[i][j]=-U[i][j+1];
                    T[i][j]=T[i][j+1];
                    break;
                case 130: //bw
                    U[i-1][j]=0;
                    V[i][j-1]=-V[i-1][j-1];
                    V[i][j]=-V[i-1][j];
                    T[i][j]=T[i-1][j];
                    break;
                case 66: //bs
                    V[i][j-1]=0;
                    U[i-1][j]=-U[i-1][j-1];
                    U[i][j]=-U[i][j-1];
                    T[i][j]=T[i][j-1];
                    break;
                case 258: //bo
                    U[i][j]=0;
                    V[i][j]=-V[i+1][j];
                    V[i][j-1]=-V[i+1][j-1];
                    T[i][j]=T[i+1][j];
                    break;
                    
                    //corner cells
                case 290: //BN, BO
                    U[i][j]=0;
                    V[i][j]=0;
                    U[i-1][j]=-U[i-1][j+1];
                    V[i][j-1]=-V[i+1][j-1];
                    T[i][j]=0.5*(T[i+1][j]+T[i][j+1]);
                    break;
                    
                case 162: //bn, bw
                    U[i-1][j]=0;
                    V[i][j]=0;
                    U[i][j]=-U[i][j+1];
                    V[i][j-1]=-V[i-1][j-1];
                    T[i][j]=0.5*(T[i-1][j]+T[i][j+1]);
                    break;
                    
                case 194: //bw, bs
                    U[i-1][j]=0;
                    V[i][j-1]=0;
                    U[i][j]=-U[i][j-1];
                    V[i][j]=-V[i-1][j];
                    T[i][j]=0.5*(T[i-1][j]+T[i][j-1]);
                    break;
                    
                case 322: //bs, bo
                    U[i][j]=0;
                    V[i][j-1]=0;
                    U[i-1][j]=-U[i-1][j-1];
                    V[i][j]=-V[i+1][j];
                    T[i][j]=0.5*(T[i+1][j]+T[i][j-1]);
                    break;
                    
                    //      default: //need to change this to show error
                    //        break; //this is for fluids
            }
            j++;
        }
        i++;
    }
}



/*

 boundary_val
 No-slip 258, 130, 66, 34, 290, 322, 194, 162
 Free-slip 260, 132, 68, 36, 292, 324, 196, 164
 outflow 264, 136, 72, 40, 296, 328, 200, 168
 inflow 272, 144, 80, 48, 304, 336, 208, 176

 BN 34
 BW 130
 BS 66
 BO 258
 
 BNO 290
 BNW 162
 BWS 194
 BSO 322
 */
