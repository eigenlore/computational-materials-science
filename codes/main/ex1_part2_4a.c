
/*******************************************************************************
 *
 * File ex1_part2_4a.c
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
    FILE *fd1, *fd2, *fd3, *fd4;

    sprintf(file_name, "../data/input_files/fcc100a%d.dat", N);
    load_data(file_name);

    thermalization();

    sprintf(file_name, "../data/ex1_part2/4a/energy_temperature.dat");
    fd1 = fopen(file_name, "w");
    sprintf(file_name, "../data/ex1_part2/4a/x_positions.dat");
    fd2 = fopen(file_name, "w");
    sprintf(file_name, "../data/ex1_part2/4a/y_positions.dat");
    fd3 = fopen(file_name, "w");
    sprintf(file_name, "../data/ex1_part2/4a/z_positions.dat");
    fd4 = fopen(file_name, "w");

    for (i = 0; i * DT < TOT_TIME; i++)
    {
        fprintf(fd1, "%.15e %.15e %.15e\n", i * DT, eval_K() + eval_U(), eval_temperature());

        for (j = 0; j < N; j++)
        {
            fprintf(fd2, "%.15e ", xx[j]);
            fprintf(fd3, "%.15e ", yy[j]);
            fprintf(fd4, "%.15e ", zz[j]);
        }
        fprintf(fd2, "\n");
        fprintf(fd3, "\n");
        fprintf(fd4, "\n");

        verlet_evolution();
    }

    fclose(fd1);
    fclose(fd2);
    fclose(fd3);
    fclose(fd4);



    free_all();

    return 0;
}
