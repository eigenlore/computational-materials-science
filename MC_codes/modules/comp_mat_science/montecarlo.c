
/*******************************************************************************
 *
 * Library montecarlo.c
 *
 * The externally accessible functions are:
 *
 *
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "global.h"
#include "random.h"
#include "montecarlo.h"

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

void init_configuration()
{
    int i, j, x, y;
    double r;

    eval_list_nbrs();

    rlxd_init(1, seed);

    for (i = 0; i < LX; i++)
        for (j = 0; j < LY; j++)
            occupation_matrix[i][j] = 0;

    i = 0;

    while (i < N)
    {
        ranlxd(&r, 1);
        x = (int)(r * LX);
        ranlxd(&r, 1);
        y = (int)(r * LY);

        assert(x < LX);
        assert(y < LY);

        if (occupation_matrix[x][y] == 0)
        {
            occupation_matrix[x][y] = 1;
            atom_position[i][0] = x;
            atom_position[i][1] = y;
            i++;
        }
    }
}

void eval_list_nbrs()
{
    int i, j;

    for (i = 0; i < LX; i++)
    {
        for (j = 0; j < LY; j++)
        {
            left_nbrs[i][j][0] = (i - 1 + LX) % LX;
            left_nbrs[i][j][1] = j;
            right_nbrs[i][j][0] = (i + 1 + LX) % LX;
            right_nbrs[i][j][1] = j;
            up_nbrs[i][j][0] = i;
            up_nbrs[i][j][1] = (j + 1 + LY) % LY;
            down_nbrs[i][j][0] = i;
            down_nbrs[i][j][1] = (j - 1 + LY) % LY;
        }
    }
}

void print_configuration(char file_name[])
{
    int i;
    FILE *fd;

    fd = fopen(file_name, "w");

    for (i = 0; i < N; i++)
        fprintf(fd, "%d %d\n", atom_position[i][0], atom_position[i][1]);

    fclose(fd);
}

int number_of_nbrs(int i, int j)
{
    return occupation_matrix[left_nbrs[i][j][0]][left_nbrs[i][j][1]] +
           occupation_matrix[right_nbrs[i][j][0]][right_nbrs[i][j][1]] +
           occupation_matrix[up_nbrs[i][j][0]][up_nbrs[i][j][1]] +
           occupation_matrix[down_nbrs[i][j][0]][down_nbrs[i][j][1]];
}

double eval_E()
{
    int i;
    double E;

    E = 0;

    for (i = 0; i < N; i++)
        E += 0.5 * J1 * number_of_nbrs(atom_position[i][0], atom_position[i][1]);

    return E;
}

double mean_number_of_nbrs()
{
    int i;
    double mean;

    mean = 0;

    for (i = 0; i < N; i++)
        mean += number_of_nbrs(atom_position[i][0], atom_position[i][1]);

    mean /= (double)N;

    return mean;
}

void sweep()
{
    int x, y, atom, x_old, y_old;
    double r, E_old;

    /*Select a random atom*/
    ranlxd(&r, 1);
    atom = (int)(r * N);
    assert(atom < N);

    x_old = atom_position[atom][0];
    y_old = atom_position[atom][1];
    E_old = eval_E();

    /*Move the atom to a new random empty slot*/
    while (1) /*repeat until it finds an empty slot*/
    {
        ranlxd(&r, 1);
        x = (int)(r * LX);
        ranlxd(&r, 1);
        y = (int)(r * LY);

        if (occupation_matrix[x][y] == 0)
        {
            occupation_matrix[x][y] = 1;
            occupation_matrix[x_old][y_old] = 0;
            atom_position[atom][0] = x;
            atom_position[atom][1] = y;
            break;
        }
    }

    if (eval_E() - E_old < 1e-8) /*E_new <= E_old*/
    {
        /*accepts the new configuration*/
        return;
    }
    else /*E_new > E_old*/
    {
        if (T == 0)
        {
            /*rejects the new configuration*/
            occupation_matrix[x][y] = 0;
            occupation_matrix[x_old][y_old] = 1;
            atom_position[atom][0] = x_old;
            atom_position[atom][1] = y_old;
            return;
        }
        else /*T not 0*/
        {
            ranlxd(&r, 1);
            if (r < exp((E_old - eval_E()) / (KB * T)))
            {
                /*accepts the new configuration*/
                return;
            }
            else
            {
                /*rejects the new configuration*/
                occupation_matrix[x][y] = 0;
                occupation_matrix[x_old][y_old] = 1;
                atom_position[atom][0] = x_old;
                atom_position[atom][1] = y_old;
                return;
            }
        }
    }
}

void thermalization(char file_name[])
{
    int i;
    FILE *fd;

    fd = fopen(file_name, "w");

    eval_list_nbrs();
    init_configuration();

    for (i = 0; i < N_TERM; i++)
    {
        fprintf(fd, "%.15e\n", eval_E());
        sweep();
    }

    fclose(fd);
}