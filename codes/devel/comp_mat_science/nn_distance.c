
/*******************************************************************************
 *
 * File nn_distance.c
 *
 * Test to check the functions load_data and eval_nn_distance. This code 
 * evaluates the nn distance of fcc100a256. 
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

    return 0;
}
