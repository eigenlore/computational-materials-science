
/*******************************************************************************
 *
 * File print_potential.c
 *
 * Print on file the Lennard Jones potential.
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
    print_potential();

    return 0;
}
