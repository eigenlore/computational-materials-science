#ifndef MONTECARLO_H
#define MONTECARLO_H

double powerd(double x, int y);
void init_configuration();
void eval_list_nbrs();
void print_configuration(char file_name[]);
double eval_E();
double mean_number_of_nbrs();
void sweep();
void thermalization(char file_name[]);

#ifdef DIM2
int number_of_nbrs(int i, int j);
#endif /*DIM2*/

#ifdef DIM3
int number_of_nbrs(int i, int j, int k);
int count_first_layer();
void init_configuration_first_layer();
void thermalization_first_layer(char file_name[]);
#endif /*DIM3*/

#endif /*MONTECARLO_H*/
