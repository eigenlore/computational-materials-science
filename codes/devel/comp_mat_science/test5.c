
/*******************************************************************************
 *
 * File test4.c
 *
 * Test to check the Verlet evolution
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
    generate_inital_v(T_INIT);
    eval_forces();

    for(i=0; i<NSTEPS; i++)
    {
        printf("The energy is %f\t", eval_K()+eval_U());
        printf("The temperature is %f\n", eval_temperature());
        verlet_evolution();
    }

    

    return 0;
}
