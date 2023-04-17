#include "utils.h"

// Initiate random matrix or a certain IC
void init_halo_mat(int *mat, int halo_sz) {
  for (int i = 0; i < halo_sz; i++) {
    for (int j = 0; j < halo_sz; j++) {
      mat[j+i*halo_sz] = rand()%2;
      // mat[j+i*halo_sz] = j+i*halo_sz;

      // Two gliders for testing periodic BC:s
      // TWO_GLIDERS
    }
  }
}

// Print a matrix of a certain dimension
void print_mat(int *mat, int mat_dim) {
  for (int i = 0; i < mat_dim+1; ++i) {
    printf("--");
  }
  printf("\n");
  for (int i = 0; i < mat_dim; i++) {
    printf("|");
    for (int j = 0; j < mat_dim; j++) {
      if (mat[j+i*mat_dim] == 1) {
        printf("\u2B23 ");
      }
      else
        printf("  ");
      // printf("%d\t", mat[j+i*mat_dim]);
    }
    printf("|\n");
  }
  for (int i = 0; i < mat_dim+1; ++i) {
    printf("--");
  }
  printf("\n");
}

// Game of Life stencil and rules
void game_of_life(int halo_sz, int *state, int *next_state) {
  int nbr_cells, state_pos = 0;
  for (int i = 0; i < halo_sz-2; i++) {
    for (int j = halo_sz+1; j < 2*halo_sz-1; j++) {
      // 2D 8-point stencil in row-major memory
      // where [i+j*halo_sz] is the center cells
      // inside the halo region.
      nbr_cells = state[j+i*halo_sz-halo_sz-1]
        + state[j+i*halo_sz-halo_sz]
        + state[j+i*halo_sz-halo_sz+1]
        + state[j+i*halo_sz-1]
        + state[j+i*halo_sz+1]
        + state[j+i*halo_sz+halo_sz-1]
        + state[j+i*halo_sz+halo_sz]
        + state[j+i*halo_sz+halo_sz+1];

      if (nbr_cells == 3 || state[j+i*halo_sz] + nbr_cells == 3) {
        next_state[j+i*halo_sz] = 1;
      }
      else
        next_state[j+i*halo_sz] = 0;
      state_pos++;
    }
  }
}
