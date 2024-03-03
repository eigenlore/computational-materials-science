
/*******************************************************************************
 *
 * File global.h
 *
 * Global parameters and arrays
 *
 *
 *
 *
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define KB 0.00008618460742911316 /*eV/K*/
#define LX 25
#define LY 25
#define N 25
#define J1 -0.2 /*eV*/
#define N_SWEEP 500000
#define T 2000
#define N_TERM 200

#ifdef MAIN_PROGRAM
#define EXTERN
#else
#define EXTERN extern
#endif /*MAIN_PROGRAM*/

EXTERN short occupation_matrix[LX][LY];
EXTERN int left_nbrs[LX][LY][2];
EXTERN int right_nbrs[LX][LY][2];
EXTERN int up_nbrs[LX][LY][2];
EXTERN int down_nbrs[LX][LY][2];
EXTERN int atom_position[N][2];
EXTERN int seed;

#undef EXTERN

#endif /*GLOBAL_H*/
