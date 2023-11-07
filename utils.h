#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>
#include <SDL2/SDL.h>

// Args and errors
#define NUMBER_OF_ARGS 3
#define INIT_ERR 1
#define PERFECT_SQUARE_ERR 2

// For the cartesian topology
#define NEAREST_NBRS 1
#define NBRS 4
#define ROOT 0
#define TAG 0
#define COL_DIM 1
#define ROW_DIM 0

// Directions in the matrix are UP, DOWN, LEFT, RIGHT 
// and directions along comms NORTH, SOUTH, WEST, EAST
#define UP_ROW_POS halo_sz+1
#define DOWN_ROW_POS block_sz*halo_sz+1
#define UP_HALO_POS 1
#define DOWN_HALO_POS block_sz*halo_sz+halo_sz+1
#define RIGHT_COL_POS halo_sz-2
#define LEFT_HALO_POS 0
#define RIGHT_HALO_POS halo_sz-1
#define LEFT_COL_POS 1

// Alternating offset for putting BC in this or next timestep
#define TIMESTEP_OFFSET halo_sz*halo_sz*((i+1)%2)

// Clears the terminal and prints the computed state from ROOT
#define CLEAR_AND_PRINT_STATE usleep(40000); if (cart_rank == ROOT) \
{printf("\e[1;1H\e[2J\n"); print_mat(this_state, halo_sz);}

// IC locations used below
#define OFS1 (mat_dim*mat_dim)/2+mat_dim/3
#define OFS2 (mat_dim*mat_dim)/4-mat_dim/3

// IC with two gliders in opposite direction to test toroidal BC
#define TWO_GLIDERS \
      mat[OFS1] = 1; \
      mat[OFS1-1] = 1; \
      mat[OFS1+mat_dim] = 1; \
      mat[OFS1+mat_dim-2] = 1; \
      mat[OFS1+2*mat_dim] = 1; \
      mat[OFS2] = 1; \
      mat[OFS2+1] = 1; \
      mat[OFS2+mat_dim] = 1; \
      mat[OFS2+2*mat_dim] = 1; \
      mat[OFS2+mat_dim+2] = 1;

void init_halo_mat(int *mat, int mat_dim);
void game_of_life(int halo_sz, int *state, int *next_state);
void print_mat(int *mat, int halo_sz);

#endif /* UTILS_H_ */ 
