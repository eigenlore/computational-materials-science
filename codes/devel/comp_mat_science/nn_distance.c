
/*******************************************************************************
 *
 * File test1.c
 *
 * Test to check the functions load_data, eval_nn_distance and eval_U. This code 
 * evaluates the nn distance and potential energy of fcc100a256. 
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

int main(int argc, char *argv[])
{
    char input_file_name[100];

    sprintf(input_file_name, "../../data/input_files/fcc100a%d.dat", N);
    load_data(input_file_name);
    printf("The nn distance is %f\n", eval_nn_distance());
    printf("The total potential is %f\n", eval_U());

    return 0;
}
