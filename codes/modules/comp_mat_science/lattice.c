
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

double eval_dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

double eval_nn_distance(double *x, double *y, double *z)
{
    double nn_distance, temp_dist;
    int i, j;

    nn_distance = eval_dist(x[0], y[0], z[0], x[1], y[1], z[1]); /*first guess*/

    for (i = 0; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            temp_dist = eval_dist(x[i], y[i], z[i], x[j], y[j], z[j]);
            if (temp_dist < nn_distance)
                nn_distance = temp_dist;
        }
    }

    return nn_distance;
}
