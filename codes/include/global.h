
/*******************************************************************************
 *
 * File global.h
 *
 * Global parameters and arrays
 *
 * N number of atoms
 * EPS, SIGMA parameters of Lennard Jones potential (eV and Angstrom)
 * RC cutoff radius for Lennard Jones
 * x, y, z atoms positions in the lattice
 *
 *******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define N 256
#define EPS 0.345
#define SIGMA 2.644
#define RC 4.5

#ifdef MAIN_PROGRAM
#define EXTERN
#else
#define EXTERN extern
#endif /*MAIN_PROGRAM*/

EXTERN double xx[N];
EXTERN double yy[N];
EXTERN double zz[N];

#undef EXTERN

#endif /*GLOBAL_H*/
