
/*******************************************************************************
 *
 * File ex1_part1_1d.c
 *
 * Print on file energy and temperature at each time step.
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/
#define MAIN_PROGRAM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "lattice.h"
#include "random.h"
#include <assert.h>

int main(int argc, char *argv[])
{
    int i;
    char file_name[100];
    FILE *fd;

    sprintf(file_name, "../data/input_files/fcc100a%d.dat", N);
    load_data(file_name);

    sprintf(file_name, "../data/ex1_part1/1d/force_and_U.dat");
    steepest_descent(file_name);

    sprintf(file_name, "../data/ex1_part1/1d/energy_temperature.dat");
    fd = fopen(file_name, "w");

    sprintf(file_name, "../data/ex1_part1/1d/therm_energy_temperature.dat");
    thermalization(file_name);

    for (i = 0; i * DT < TOT_TIME; i++)
    {
        fprintf(fd, "%.15e %.15e %.15e\n", i * DT, eval_K() + eval_U(), eval_temperature());
        verlet_evolution();
    }

    fclose(fd);

    free_all();

    return 0;
}
