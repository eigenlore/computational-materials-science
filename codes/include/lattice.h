#ifndef LATTICE_H
#define LATTICE_H

void load_data(double *x, double *y, double *z, char file_name[]);
double eval_dist(double x1, double y1, double z1, double x2, double y2, double z2);
double eval_nn_distance(double *x, double *y, double *z);



#endif /*LATTICE_H*/
