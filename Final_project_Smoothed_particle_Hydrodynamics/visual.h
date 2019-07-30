    #ifndef __VISUAL_H__
    #define __VISUAL_H__

    void write_vtkFile(const char *szProblem,
        int timeStepNumber,
        double xlength,
    double ylength,
    double x_origin,
    double y_origin,
    int num_inner_particles,
    int num_boundary_particles,
        double dx,
        double dy,
    struct particle* domainparticles,
    int temp);


    void write_vtkHeader( FILE *fp, int num_inner_particles, int num_boundary_particles,
                      double dx, double dy,double x_length, double y_length, int temp);


    void write_vtkPointCoordinates( FILE *fp, int num_inner_particles, int num_boundary_particles,
                      double dx, double dy, struct particle* domainparticles, int temp);          




    #endif
