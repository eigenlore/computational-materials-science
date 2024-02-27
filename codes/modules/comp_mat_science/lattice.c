
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
 *      potential. It sums only on neighbors.
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

double powerd(double x, int y)
{
    double temp;
    if (y == 0)
        return 1;
    temp = powerd(x, y / 2);
    if ((y % 2) == 0)
    {
        return temp * temp;
    }
    else
    {
        if (y > 0)
            return x * temp * temp;
        else
            return (temp * temp) / x;
    }
}

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

static double eval_dist1D(double a, double b)
{
    double res;

    res = a - b;
    if (PBC)
        res -= SIZE * floor(res / SIZE + 0.5);

    return res;
}

static double eval_dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt(eval_dist1D(x1, x2) * eval_dist1D(x1, x2) + eval_dist1D(y1, y2) * eval_dist1D(y1, y2) + eval_dist1D(z1, z2) * eval_dist1D(z1, z2));
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
        return 4 * EPS * (powerd(SIGMA / r, 12) - powerd(SIGMA / r, 6));
    else
        return 0;
}

double eval_U()
{
    int i, j, k;
    double U;

    U = 0;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < number_nbrs[i]; j++)
        {
            k = which_nbrs[i][j];
            U += lennard_jones(eval_dist(xx[i], yy[i], zz[i], xx[k], yy[k], zz[k]));
        }
    }
    U /= 2;

    return U;
}

void eval_nbrs()
{
    int i, j, count, *temp;

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

void generate_inital_v()
{
    int i;
    double C, r[3 * N], v_tot_x, v_tot_y, v_tot_z, T_temp;

    rlxd_init(1, 3122000); /*seed*/
    ranlxd(r, 3 * N);
    C = sqrt(3 * KB * T_INIT / M);
    v_tot_x = 0;
    v_tot_y = 0;
    v_tot_z = 0;

    /*Generate inital velocities*/
    for (i = 0; i < N; i++)
    {
        vxx[i] = 2 * C * (r[i] - 0.5);
        vyy[i] = 2 * C * (r[i + N] - 0.5);
        vzz[i] = 2 * C * (r[i + 2 * N] - 0.5);
        v_tot_x += vxx[i];
        v_tot_y += vyy[i];
        v_tot_z += vzz[i];
    }

    /*Adjust to be to have stationary center of mass*/
    for (i = 0; i < N; i++)
    {
        vxx[i] -= v_tot_x / (double)N;
        vyy[i] -= v_tot_y / (double)N;
        vzz[i] -= v_tot_z / (double)N;
    }

    T_temp = eval_temperature();

    /*Riscale to have temperature T instead of T_temp*/
    for (i = 0; i < N; i++)
    {
        vxx[i] *= sqrt(T_INIT / T_temp);
        vyy[i] *= sqrt(T_INIT / T_temp);
        vzz[i] *= sqrt(T_INIT / T_temp);
    }
}

double eval_K()
{
    int i;
    double K;
    K = 0;

    for (i = 0; i < N; i++)
        K += vxx[i] * vxx[i] + vyy[i] * vyy[i] + vzz[i] * vzz[i];

    K *= M / 2;
    return K;
}

double eval_temperature()
{
    return 2 * eval_K() / (3 * N * KB);
}

void eval_forces()
{
    int i, j, k;
    double r;

    for (i = 0; i < N; i++)
    {
        Fxx[i] = 0;
        Fyy[i] = 0;
        Fzz[i] = 0;
        for (j = 0; j < number_nbrs[i]; j++)
        {
            k = which_nbrs[i][j];
            r = eval_dist(xx[i], yy[i], zz[i], xx[k], yy[k], zz[k]);
            Fxx[i] += 24 * EPS * powerd(SIGMA, 6) * powerd(r, -8) * eval_dist1D(xx[i], xx[k]) * (2 * powerd(SIGMA, 6) * powerd(r, -6) - 1);
            Fyy[i] += 24 * EPS * powerd(SIGMA, 6) * powerd(r, -8) * eval_dist1D(yy[i], yy[k]) * (2 * powerd(SIGMA, 6) * powerd(r, -6) - 1);
            Fzz[i] += 24 * EPS * powerd(SIGMA, 6) * powerd(r, -8) * eval_dist1D(zz[i], zz[k]) * (2 * powerd(SIGMA, 6) * powerd(r, -6) - 1);
        }
    }
}

void verlet_evolution()
{
    int i;
    double old_Fx[N], old_Fy[N], old_Fz[N];

    for (i = 0; i < N; i++)
    {
        old_Fx[i] = Fxx[i];
        old_Fy[i] = Fyy[i];
        old_Fz[i] = Fzz[i];

        xx[i] += vxx[i] * DT + Fxx[i] * DT * DT / (2 * M);
        yy[i] += vyy[i] * DT + Fyy[i] * DT * DT / (2 * M);
        zz[i] += vzz[i] * DT + Fzz[i] * DT * DT / (2 * M);
    }

    eval_nbrs();
    eval_forces();

    for (i = 0; i < N; i++)
    {
        vxx[i] += (Fxx[i] + old_Fx[i]) * DT / (2 * M);
        vyy[i] += (Fyy[i] + old_Fy[i]) * DT / (2 * M);
        vzz[i] += (Fzz[i] + old_Fz[i]) * DT / (2 * M);
    }
}

void thermalization()
{
    int i;

    eval_nbrs();
    eval_forces();
    generate_inital_v();

    for (i = 0; i * DT < TERM_TIME; i++)
        verlet_evolution();
}