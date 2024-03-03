
/*******************************************************************************
 *
 * File general_test.c
 *
 * General tests function.
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/
#define MAIN_PROGRAM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "montecarlo.h"
#include "random.h"
#include <assert.h>
#include <time.h>

int main(int argc, char *argv[])
{
    seed = time(NULL);

    printf("%f\n", powerd(2, 31) - 1);
    printf("%d\n", seed);

    return 0;
}