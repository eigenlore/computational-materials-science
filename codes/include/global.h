
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
#define RP 4.2                    /*A*/
#define KB 0.00008618460742911316 /*eV/K*/
#define M 11.205e-27              /*kg*/
#define DT 1e-15                  /*seconds*/
#define T_INIT 1550               /*Kelvin*/
#define TERM_TIME 3e-12           /*seconds*/
#define TOT_TIME 10e-12           /*seconds*/
#define PBC 1                     /*1 with PBC, 0 without*/
#define SIZE 16.641600            /*A*/
#define A 1.762150591975983e+08
#define B -2.840236130038733e+08
#define C 1.961490151653693e+08
#define D -7.523885828765376e+07
#define E 1.731196189189488e+07
#define F -2.389452170224014e+06
#define G 1.831785386120462e+05
#define H -6.016872195722493e+03

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
