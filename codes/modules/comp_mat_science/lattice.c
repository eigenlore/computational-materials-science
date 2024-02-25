
/*******************************************************************************
 *
 * Library lattice.c
 *
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "global.h"
#include "random.h"
#include "lattice.h"

void load_data(double *x, double *y, double *z, char file_name[])
{
    FILE *file;
    int row;

    file = fopen(file_name, "r");
    assert(file != NULL);

    for (row = 0; row < N; row++)
        assert(fscanf(file, "%lf %lf %lf", x + row, y + row, z + row) == 3);
}




double eval_nn_distance(double *x, double *y, double *z)
{
}
