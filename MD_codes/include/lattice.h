#ifndef LATTICE_H
#define LATTICE_H

double powerd(double x, int y);
void free_all();
void load_data(char file_name[]);
double eval_nn_distance();
double eval_U();
void eval_nbrs();
void generate_inital_v();
double eval_K();
double eval_temperature();
void eval_forces();
void verlet_evolution();
void euler_evolution();
void thermalization(char file_name[]);
void eval_coefficients();
void print_potential();
double *eval_L();
double *eval_v_cm();
double eval_max_force();
void steepest_descent(char file_name[]);

#endif /*LATTICE_H*/
