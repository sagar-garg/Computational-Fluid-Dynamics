#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
    int imax,
    int jmax,
    double **U,
    double **V, int **flag, double UI, double VI,double U_in, double V_in) ;

void spec_boundary_val (int imax, int jmax, double **U, double **V, int **flag, double UI, double VI);


//For temperature//
void boundaryvalues_t(
                    int imax,
                    int jmax,
                    double **U,
                    double **V,double **T, int **flag, double UI, double VI,double U_in, double V_in, double Te,double Tw, double Tn, double Ts);

void spec_boundary_val_t(int imax, int jmax, double **U, double **V, int **flag, double **T);

#endif

