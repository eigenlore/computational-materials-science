
/*******************************************************************************
 *
 * File test4.c
 *
 * Test to check the functions eval_forces
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
    int i;
    char input_file_name[100];

    sprintf(input_file_name, "../../data/input_files/fcc100a%d.dat", N);
    load_data(input_file_name);
    eval_nbrs();
    eval_forces();
    for (i = 0; i < N; i++)
    {
        printf("%d %d %f %f %f\n", i + 1, number_nbrs[i], Fxx[i], Fyy[i], Fzz[i]);
    }

    return 0;
}
