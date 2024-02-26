#ifndef LATTICE_H
#define LATTICE_H

double powerd(double x, int y);
void load_data(char file_name[]);
double eval_nn_distance();
double eval_U();
void eval_nbrs();
void generate_inital_v(double T);
double eval_K();
double eval_temperature();
void eval_forces();
void verlet_evolution();

#endif /*LATTICE_H*/
