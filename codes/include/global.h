
/*******************************************************************************
 *
 * File global.h
 *
 * Global parameters and arrays
 *
 * N number of atoms
 * x, y, z atoms positions in the lattice
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

EXTERN double xx[N];
EXTERN double yy[N];
EXTERN double zz[N];

#undef EXTERN



#endif /*GLOBAL_H*/
