
/*******************************************************************************
 *
 * Library lattice.c
 *
 * The externally accessible functions are
 *
 *  void load_data(double *x, double *y, double *z, char file_name[])
 *      Load the data from the file file_name into x, y and z. 
 *
 *  double eval_nn_distance(double *x, double *y, double *z)
 *      Evaluates the nn distance of a lattice with positions of the atomes 
 *      given in x, y and z
 *
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

void load_data(char file_name[])
{
    FILE *file;
    int row;

    file = fopen(file_name, "r");
    assert(file != NULL);

    for (row = 0; row < N; row++)
        assert(fscanf(file, "%lf %lf %lf", xx + row, yy + row, zz + row) == 3);
}

double eval_dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

double eval_nn_distance()
{
    double nn_distance, temp_dist;
    int i, j;

    nn_distance = eval_dist(xx[0], yy[0], zz[0], xx[1], yy[1], zz[1]); /*first guess*/

    for (i = 0; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            temp_dist = eval_dist(xx[i], yy[i], zz[i], xx[j], yy[j], zz[j]);
            if (temp_dist < nn_distance)
                nn_distance = temp_dist;
        }
    }

    return nn_distance;
}
