
/*******************************************************************************
 *
 * File test2.c
 *
 * It writes on file how many atoms have a certain number of neighbors. 
 * It defines a class of equivalency. Two atoms are equivalent if they have
 * the same number of neighbors. 
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

#define NNFCC 12

int main(int argc, char *argv[])
{
    char file_name[100];
    int i, equivalent[NNFCC];
    FILE *file;

    sprintf(file_name, "../../data/input_files/fcc100a%d.dat", N);
    load_data(file_name);

    eval_nbrs();
   

    for (i = 0; i < NNFCC; i++)
        equivalent[i] = 0;

    for (i = 0; i < N; i++)
    {
        assert(number_nbrs[i] <= NNFCC);
        equivalent[number_nbrs[i]-1]++;
        if(number_nbrs[i]==3)
            printf("%d\n", i+1);
    }

    sprintf(file_name, "../../data/test_files/equivalent.dat");
    file = fopen(file_name, "w");
    for (i = 0; i < NNFCC; i++)
        fprintf(file, "%d\t%d\n", i+1, equivalent[i]);
    fclose(file);

   

    return 0;
}
