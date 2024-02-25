
/*******************************************************************************
 *
 * File nn_distance.c
 *
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "lattice.h"
#include "random.h"

int main(int argc, char *argv[])
{
    char input_file_name[100];
    double *x, *y, *z;

    sprintf(input_file_name, "../../data/input_files/fcc100a%d.dat", N);
    x = (double *)malloc(N * sizeof(double));
    y = (double *)malloc(N * sizeof(double));
    z = (double *)malloc(N * sizeof(double));
    load_data(x, y, z, input_file_name);

    printf("The nn distance is %f\n", eval_nn_distance(x, y, z));

    return 0;
}
