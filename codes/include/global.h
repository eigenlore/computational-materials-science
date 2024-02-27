
/*******************************************************************************
 *
 * File global.h
 *
 * Global parameters and arrays
 *
 * N number of atoms
 * EPS, SIGMA parameters of Lennard Jones potential
 * RC cutoff radius for Lennard Jones
 * x, y, z atoms positions in the lattice
 *
 *
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define N 256                     /*atoms*/
#define EPS 0.345                 /*eV*/
#define SIGMA 2.644               /*A*/
#define RC 4.5                    /*A*/
#define KB 0.00008618460742911316 /*eV/K*/
#define M 11.205e-27              /*kg*/
#define DT 8e-15                  /*seconds*/
#define T_INIT 15                 /*Kelvin*/
#define TERM_TIME 3e-12           /*seconds*/
#define TOT_TIME 10e-12           /*seconds*/
#define PBC 0                     /*1 with PBC, 0 without*/
#define SIZE 16.641600            /*A*/

#ifdef MAIN_PROGRAM
#define EXTERN
#else
#define EXTERN extern
#endif /*MAIN_PROGRAM*/

EXTERN double xx[N];
EXTERN double yy[N];
EXTERN double zz[N];
EXTERN int number_nbrs[N];
EXTERN int *which_nbrs[N];
EXTERN double vxx[N];
EXTERN double vyy[N];
EXTERN double vzz[N];
EXTERN double Fxx[N];
EXTERN double Fyy[N];
EXTERN double Fzz[N];

#undef EXTERN

#endif /*GLOBAL_H*/
