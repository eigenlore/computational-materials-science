
/*******************************************************************************
 *
 * File global.h
 *
 * Global parameters and arrays
 *
 * N number of atoms
 *
 *******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define N 256

#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif /*MAIN_PROGRAM*/

EXTERN double x[N];
EXTERN double y[N];
EXTERN double z[N];

#undef EXTERN



#endif /*GLOBAL_H*/
