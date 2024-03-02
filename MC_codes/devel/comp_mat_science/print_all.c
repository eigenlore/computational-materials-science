
/*******************************************************************************
 *
 * File print_all.c
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
    int i, count;
    char file_name[100];
    FILE *fd1, *fd2, *fd3;
    double *temp;
    sprintf(file_name, "../../data/input_files/fcc100a%d.dat", N);
    load_data(file_name);

    thermalization();

    sprintf(file_name, "../../data/test_files/energy_temperaturePBC.dat");
    fd1 = fopen(file_name, "w");
    sprintf(file_name, "../../data/test_files/angular_momentumPBC.dat");
    fd2 = fopen(file_name, "w");
    sprintf(file_name, "../../data/test_files/v_cmPBC.dat");
    fd3 = fopen(file_name, "w");

    count = 0;

    for (i = 0; i * DT < TOT_TIME; i++)
    {

        if (i == count)
        {
            printf("Completed %.0f%%.\r", i * DT * 100 / TOT_TIME);
            fflush(stdout);
            count = count + TOT_TIME / (DT * 100);
        }

        fprintf(fd1, "%.15e %.15e %.15e\n", i * DT, eval_K() + eval_U(), eval_temperature());

        temp = eval_L();
        fprintf(fd2, "%.15e %.15e %.15e\n", temp[0], temp[1], temp[2]);
        free(temp);

        temp = eval_v_cm();
        fprintf(fd3, "%.15e %.15e %.15e\n", temp[0], temp[1], temp[2]);
        free(temp);

        verlet_evolution();
    }

    fclose(fd1);
    fclose(fd2);
    fclose(fd3);

    free_all();

    return 0;
}
