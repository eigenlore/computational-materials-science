#ifndef LATTICE_H
#define LATTICE_H

void load_data(char file_name[]);
double eval_nn_distance();
double eval_U();
void eval_nbrs(int *number_nbrs, int **which_nbrs);

#endif /*LATTICE_H*/
