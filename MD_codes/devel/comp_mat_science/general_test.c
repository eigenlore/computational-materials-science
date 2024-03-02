#define MAIN_PROGRAM
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int *ptr;

int main()
{
	printf("%p\n", ptr);
	ptr=malloc(10*sizeof(int));
	free(ptr);
	printf("%p\n", ptr);
	ptr=malloc(10*sizeof(int));
	printf("%p\n", ptr);



	return 0;
}
