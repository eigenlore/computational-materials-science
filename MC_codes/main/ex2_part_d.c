
/*******************************************************************************
 *
 * File ex2_part_d.c
 *
 * General tests function.
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/
#define MAIN_PROGRAM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "montecarlo.h"
#include "random.h"
#include <assert.h>
#include <time.h>

int main(int argc, char *argv[])
{
    int i;
    char file_name[100];
    FILE *fd;

    if (argc == 2)
        seed = atoi(argv[1]);
    else
        seed = time(NULL);

    sprintf(file_name, "../data/ex2_part_d/seedJ0%.1fT%d.dat", J0, T);
    fd = fopen(file_name, "w");
    fprintf(fd, "%d\n", seed);
    fclose(fd);

    sprintf(file_name, "../data/ex2_part_d/thermalization_energyJ0%.1fT%d.dat", J0, T);

    thermalization(file_name);

    sprintf(file_name, "../data/ex2_part_d/energy_and_nbrsJ0%.1fT%d.dat", J0, T);
    fd = fopen(file_name, "w");

    for (i = 0; i < N_SWEEP; i++)
    {
        fprintf(fd, "%.15e %.15e %d\n", eval_E(), mean_number_of_nbrs(), count_first_layer());
        sweep();
    }

    fclose(fd);

    return 0;
}