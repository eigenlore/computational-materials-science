
/*******************************************************************************
 *
 * File test3.c
 *
 * Test to check the functions
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

    sprintf(input_file_name, "../../data/input_files/fcc100a%d.dat", N);
    load_data(input_file_name);

    generate_inital_v(300);
    printf("The first atom velocities are %f %f %f\n", vxx[0], vyy[0], vzz[0]);
    printf("The initial temperature is %f K\n", eval_temperature());

    return 0;
}
