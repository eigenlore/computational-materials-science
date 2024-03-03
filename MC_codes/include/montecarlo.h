#ifndef MONTECARLO_H
#define MONTECARLO_H

double powerd(double x, int y);
void init_configuration();
void eval_list_nbrs();
void print_configuration(char file_name[]);
int number_of_nbrs(int i, int j);
double eval_E();
double mean_number_of_nbrs();
void sweep();
void thermalization(char file_name[]);

#endif /*MONTECARLO_H*/
