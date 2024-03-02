
/*******************************************************************************
 *
 * File test1.c
 *
 * Test to check the functions load_data, eval_nn_distance, eval_nbrs and eval_U.
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
    char input_file_name[100];
    int i;

    sprintf(input_file_name, "../../data/input_files/fcc100a%d.dat", N);
    load_data(input_file_name);
    eval_nbrs();

    printf("The nn distance is %f\n", eval_nn_distance());
    printf("The total potential is %f eV\n", eval_U());
    printf("The energy per atom is %f eV\n", eval_U() / (double)N);
    printf("The neighbors of the third atom are:\n");
    for (i = 0; i < number_nbrs[2]; i++)
        printf("%d\n", which_nbrs[2][i] + 1); /*+1 to compare with Matlab whose indexes start from 1*/


    return 0;
}
