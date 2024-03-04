
/*******************************************************************************
 *
 * File global.h
 *
 * Global parameters and arrays
 *
 * KB Boltzmann constant
 * LX grid lenght along x
 * LY grid lenght along y
 * LZ grid lenght along z
 * N number of atoms
 * J0 energy between atoms and substrate
 * J1 energy between first neighbours atom
 * N_SWEEP number of sweeps performed after thermalization
 * T temperature
 * N_TERM number of sweeps to reach thermalization
 * DIM3 or DIM2 according to dimension
 *
 * occupation_matrix: 0 if the cell is empty, 1 otherwise
 * left_nbrs: coordinate of left neighbours
 * right_nbrs: coordinate of right neighbours
 * up_nbrs: coordinate of up neighbours
 * down_nbrs: coordinate of down neighbours
 * top_nbrs: coordinate of top neighbours
 * bottom_nbrs: coordinate of bottom neighbours
 * atom_position: coordinate of positions of all atoms
 *
 * Author: Lorenzo Tasca
 *
 *******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define KB 0.00008618460742911316 /*eV/K*/
#define LX 25
#define LY 25
#define LZ 10
#define N 27
#define J0 -0.4 /*eV*/
#define J1 -0.2 /*eV*/
#define N_SWEEP 500000
#define T 500 /*K*/
#define N_TERM 200
#define DIM3

#ifdef MAIN_PROGRAM
#define EXTERN
#else
#define EXTERN extern
#endif /*MAIN_PROGRAM*/

EXTERN int seed;

#ifdef DIM2
EXTERN short occupation_matrix[LX][LY];
EXTERN int left_nbrs[LX][LY][2];
EXTERN int right_nbrs[LX][LY][2];
EXTERN int up_nbrs[LX][LY][2];
EXTERN int down_nbrs[LX][LY][2];
EXTERN int atom_position[N][2];
#endif /*DIM2*/

#ifdef DIM3
EXTERN short occupation_matrix[LX][LY][LZ];
EXTERN int left_nbrs[LX][LY][LZ][3];
EXTERN int right_nbrs[LX][LY][LZ][3];
EXTERN int up_nbrs[LX][LY][LZ][3];
EXTERN int down_nbrs[LX][LY][LZ][3];
EXTERN int top_nbrs[LX][LY][LZ][3];
EXTERN int bottom_nbrs[LX][LY][LZ][3];
EXTERN int atom_position[N][3];
#endif /*DIM3*/

#undef EXTERN

#endif /*GLOBAL_H*/
