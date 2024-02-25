#define MAIN_PROGRAM
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <lattice.h>
#include <global.h>

int main ()
{
	int *p;
	p = (int*) malloc(N*sizeof(int));
	assert(p==NULL);
	return 0;
}
