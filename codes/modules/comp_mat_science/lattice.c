
/*******************************************************************************
 *
 * Library lattice.c
 *
 * The externally accessible functions are
 *
 *  void load_data(char file_name[])
 *      Load the data from the file file_name into global xx, yy and zz.
 *
 *  double eval_nn_distance()
 *      Evaluates the nn distance of a lattice
 *
 *  double eval_U()
 *      Evaluates the potential energy of the lattice using Lennard Jones
 *      potential
 *
 *  void eval_nbrs(int *number_nbrs, int **number_nbrs)
 *      Evaluates the number of neighbors of the atom i number_nbrs[i] and
 *      the list of their indexes which_nbrs[i]. which_nbrs[i] has lenght
 *      number_nbrs[i]. It is care of the user to allocate number_nbrs and 
 *      number_nbrs of the correct dimension and to eventually free them. 
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

    fclose(file);
}

static double eval_dist(double x1, double y1, double z1, double x2, double y2, double z2)
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

static double lennard_jones(double r)
{
    if (r < RC)
        return 4 * EPS * (pow(SIGMA / r, 12) - pow(SIGMA / r, 6));
    else
        return 0;
}

double eval_U()
{
    int i, j, k, *number_nbrs, **which_nbrs;
    double U;

    U = 0;
    number_nbrs = (int *)malloc(N * sizeof(int));
    which_nbrs = (int **)malloc(N * sizeof(int *));
    eval_nbrs(number_nbrs, which_nbrs);
    assert(number_nbrs != NULL && which_nbrs != NULL);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < number_nbrs[i]; j++)
        {
            k = which_nbrs[i][j];
            U += lennard_jones(eval_dist(xx[i], yy[i], zz[i], xx[k], yy[k], zz[k]));
        }
    }
    U/=2;

    for (i = 0; i < N; i++)
        free(which_nbrs[i]);
    free(which_nbrs);
    free(number_nbrs);

    return U;
}

void eval_nbrs(int *number_nbrs, int **which_nbrs)
{
    int i, j, count, *temp;

    assert(number_nbrs != NULL && which_nbrs != NULL); /*pointers passed must be non empty*/

    temp = (int *)malloc(N * sizeof(int));

    for (i = 0; i < N; i++)
    {
        number_nbrs[i] = 0;
        count = 0;

        /*I calculate the number of neighbors and save their indexes in temp*/
        for (j = 0; j < N; j++)
            if (eval_dist(xx[i], yy[i], zz[i], xx[j], yy[j], zz[j]) < RC && j != i)
            {
                number_nbrs[i]++;
                temp[count++] = j;
            }

        /*I allocate the memory for which_nbrs[i] and copy the indexes from temp*/
        if (number_nbrs[i] != 0)
        {
            which_nbrs[i] = (int *)malloc(number_nbrs[i] * sizeof(int));
            for (j = 0; j < number_nbrs[i]; j++)
                which_nbrs[i][j] = temp[j];
        }
    }
    free(temp);
}
