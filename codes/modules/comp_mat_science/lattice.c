
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

void free_all()
{
    int i;

    for (i = 0; i < N; i++)
    {
        free(which_nbrs[i]);
        which_nbrs[i] = NULL;
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

    if (r < RP)
        return 4 * EPS * (powerd(SIGMA / r, 12) - powerd(SIGMA / r, 6));
    else
        return A + B * powerd(r, 1) + C * powerd(r, 2) + D * powerd(r, 3) + E * powerd(r, 4) + F * powerd(r, 5) + G * powerd(r, 6) + H * powerd(r, 7);
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

    free_all();

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
    double c, r[3 * N], v_tot_x, v_tot_y, v_tot_z, T_temp;

    rlxd_init(1, 3122000); /*seed*/
    ranlxd(r, 3 * N);
    c = sqrt(3 * KB * T_INIT / M);
    v_tot_x = 0;
    v_tot_y = 0;
    v_tot_z = 0;

    /*Generate inital velocities*/
    for (i = 0; i < N; i++)
    {
        vxx[i] = 2 * c * (r[i] - 0.5);
        vyy[i] = 2 * c * (r[i + N] - 0.5);
        vzz[i] = 2 * c * (r[i + 2 * N] - 0.5);
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

void euler_evolution()
{
    int i;

    for (i = 0; i < N; i++)
    {
        xx[i] += DT * vxx[i];
        yy[i] += DT * vyy[i];
        zz[i] += DT * vzz[i];
        vxx[i] += DT * Fxx[i] / M;
        vyy[i] += DT * Fyy[i] / M;
        vzz[i] += DT * Fzz[i] / M;
    }

    eval_nbrs();
    eval_forces();
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

void eval_coefficients()
{
    double a, b, c, d, e, f, g, h;

    a = (1 / (powerd(RC - RP, 7) * powerd(RP, 12))) * 4 * EPS * powerd(RC, 4) * powerd(SIGMA, 6) * (2 * powerd(RP, 6) * (-42 * powerd(RC, 3) + 182 * powerd(RC, 2) * RP - 273 * RC * powerd(RP, 2) + 143 * powerd(RP, 3)) + (455 * powerd(RC, 3) - 1729 * powerd(RC, 2) * RP + 2223 * RC * powerd(RP, 2) - 969 * powerd(RP, 3)) * powerd(SIGMA, 6));
    b = (1 / (powerd(RC - RP, 7) * powerd(RP, 13))) * 16 * EPS * powerd(RC, 3) * powerd(SIGMA, 6) * (powerd(RP, 6) * (54 * powerd(RC, 4) - 154 * powerd(RC, 3) * RP + 351 * RC * powerd(RP, 3) - 286 * powerd(RP, 4)) + (-315 * powerd(RC, 4) + 749 * powerd(RC, 3) * RP + 171 * powerd(RC, 2) * powerd(RP, 2) - 1539 * RC * powerd(RP, 3) + 969 * powerd(RP, 4)) * powerd(SIGMA, 6));
    c = (1 / (powerd(RC - RP, 7) * powerd(RP, 14))) * 12 * EPS * powerd(RC, 2) * powerd(SIGMA, 6) * (powerd(RP, 6) * (-63 * powerd(RC, 5) - 7 * powerd(RC, 4) * RP + 665 * powerd(RC, 3) * powerd(RP, 2) - 975 * powerd(RC, 2) * powerd(RP, 3) - 52 * RC * powerd(RP, 4) + 572 * powerd(RP, 5)) + 2 * (195 * powerd(RC, 5) + 91 * powerd(RC, 4) * RP - 1781 * powerd(RC, 3) * powerd(RP, 2) + 1995 * powerd(RC, 2) * powerd(RP, 3) + 399 * RC * powerd(RP, 4) - 969 * powerd(RP, 5)) * powerd(SIGMA, 6));
    d = (1 / (powerd(RC - RP, 7) * powerd(RP, 15))) * 16 * EPS * powerd(SIGMA, 6) * (RC * powerd(RP, 6) * (14 * powerd(RC, 6) + 126 * powerd(RC, 5) * RP - 420 * powerd(RC, 4) * powerd(RP, 2) - 90 * powerd(RC, 3) * powerd(RP, 3) + 1105 * powerd(RC, 2) * powerd(RP, 4) - 624 * RC * powerd(RP, 5) - 286 * powerd(RP, 6)) + RC * (-91 * powerd(RC, 6) - 819 * powerd(RC, 5) * RP + 2145 * powerd(RC, 4) * powerd(RP, 2) + 1125 * powerd(RC, 3) * powerd(RP, 3) - 5035 * powerd(RC, 2) * powerd(RP, 4) + 1881 * RC * powerd(RP, 5) + 969 * powerd(RP, 6)) * powerd(SIGMA, 6));
    e = (1 / (powerd(RC - RP, 7) * powerd(RP, 15))) * 4 * EPS * powerd(SIGMA, 6) * (2 * powerd(RP, 6) * (-112 * powerd(RC, 6) - 63 * powerd(RC, 5) * RP + 1305 * powerd(RC, 4) * powerd(RP, 2) - 1625 * powerd(RC, 3) * powerd(RP, 3) - 585 * powerd(RC, 2) * powerd(RP, 4) + 1287 * RC * powerd(RP, 5) + 143 * powerd(RP, 6)) + (1456 * powerd(RC, 6) + 1404 * powerd(RC, 5) * RP - 14580 * powerd(RC, 4) * powerd(RP, 2) + 13015 * powerd(RC, 3) * powerd(RP, 3) + 7695 * powerd(RC, 2) * powerd(RP, 4) - 8721 * RC * powerd(RP, 5) - 969 * powerd(RP, 6)) * powerd(SIGMA, 6));
    f = (1 / (powerd(RC - RP, 7) * powerd(RP, 15))) * 48 * EPS * powerd(SIGMA, 6) * (-powerd(RP, 6) * (-28 * powerd(RC, 5) + 63 * powerd(RC, 4) * RP + 65 * powerd(RC, 3) * powerd(RP, 2) - 247 * powerd(RC, 2) * powerd(RP, 3) + 117 * RC * powerd(RP, 4) + 65 * powerd(RP, 5)) + (-182 * powerd(RC, 5) + 312 * powerd(RC, 4) * RP + 475 * powerd(RC, 3) * powerd(RP, 2) - 1140 * powerd(RC, 2) * powerd(RP, 3) + 342 * RC * powerd(RP, 4) + 228 * powerd(RP, 5)) * powerd(SIGMA, 6));
    g = (1 / (powerd(RC - RP, 7) * powerd(RP, 15))) * 4 * EPS * powerd(SIGMA, 6) * (powerd(RP, 6) * (-224 * powerd(RC, 4) + 819 * powerd(RC, 3) * RP - 741 * powerd(RC, 2) * powerd(RP, 2) - 429 * RC * powerd(RP, 3) + 715 * powerd(RP, 4)) + 2 * (728 * powerd(RC, 4) - 2223 * powerd(RC, 3) * RP + 1425 * powerd(RC, 2) * powerd(RP, 2) + 1292 * RC * powerd(RP, 3) - 1292 * powerd(RP, 4)) * powerd(SIGMA, 6));
    h = (1 / (powerd(RC - RP, 7) * powerd(RP, 15))) * 16 * EPS * powerd(SIGMA, 6) * (powerd(RP, 6) * (14 * powerd(RC, 3) - 63 * powerd(RC, 2) * RP + 99 * RC * powerd(RP, 2) - 55 * powerd(RP, 3)) + (-91 * powerd(RC, 3) + 351 * powerd(RC, 2) * RP - 459 * RC * powerd(RP, 2) + 204 * powerd(RP, 3)) * powerd(SIGMA, 6));

    printf("#define A ");
    printf("%.15e\n", a);
    printf("#define B ");
    printf("%.15e\n", b);
    printf("#define C ");
    printf("%.15e\n", c);
    printf("#define D ");
    printf("%.15e\n", d);
    printf("#define E ");
    printf("%.15e\n", e);
    printf("#define F ");
    printf("%.15e\n", f);
    printf("#define G ");
    printf("%.15e\n", g);
    printf("#define H ");
    printf("%.15e\n", h);
}