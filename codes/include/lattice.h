#ifndef LATTICE_H
#define LATTICE_H

void load_data(char file_name[]);
double eval_dist(double x1, double y1, double z1, double x2, double y2, double z2);
double eval_nn_distance();
double lennard_jones(double r);
double eval_U();



#endif /*LATTICE_H*/
