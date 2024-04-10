
/*******************************************************************************
 *
 * File ex1_part1_1abc.c
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
    int i, j;
    char file_name[100];
    FILE *fd1, *fd2;

    sprintf(file_name, "../data/input_files/fcc100a%d.dat", N);
    load_data(file_name);

    thermalization();

    sprintf(file_name, "../data/ex1_part1/1abc/velocities.dat");
    fd1 = fopen(file_name, "w");

    sprintf(file_name, "../data/ex1_part1/1abc/energy_temperature.dat");
    fd2 = fopen(file_name, "w");

    for (i = 0; i * DT < TOT_TIME; i++)
    {
        if (i % 100)
            for (j = 0; j < N; j++)
                fprintf(fd1, "%.15e %.15e %.15e\n", vxx[j], vyy[j], vzz[j]);

        fprintf(fd2, "%.15e %.15e %.15e\n", i * DT, eval_K() + eval_U(), eval_temperature());
        verlet_evolution();
    }

    fclose(fd1);
    fclose(fd2);

    free_all();

    return 0;
}
