#include "utils.h"

// Main driver function
int main(int argc, char *argv[]) {

  // Communicator size and ranks
  int comm_sz, rank, cart_rank;

  // Input args
  int args[NUMBER_OF_ARGS], N, iters, visualize;

  // Memory pointers
  int *state, *this_state, *temp_state, *next_state;
  double start, stop, time, time_max;

  MPI_Init(&argc, &argv);

  // Size of the default communicator
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Random seed dependent on rank
  srand(rank * 1337);

  // Check input args
  if (rank == ROOT) {

    // Check number of args
    if (argc != NUMBER_OF_ARGS+1) {
      printf("-----------------------------------------------------------------"
             "---------\n");
      fprintf(stderr,
          "Usage: mpirun -np <NPROCS> %s <SIZE> <ITERS> <VISUALIZE (1 or 0)>\n",
          argv[0]);
      MPI_Abort(MPI_COMM_WORLD, INIT_ERR);
    }

    // Check perfect square
    if (comm_sz % (int)sqrt(comm_sz) || atoi(argv[1]) % ((int)sqrt(comm_sz))) {
      printf("-----------------------------------------------------------------"
             "---------\n");
      fprintf(stderr, "Either <NPROCS> not a perfect square or \n<PROBLEM SIZE> not "
             "divisible by sqrt<NPROCS>\n");
      MPI_Abort(MPI_COMM_WORLD, PERFECT_SQUARE_ERR);
    }

    // Iterate through args with +1 offset 
    for (int i = 0; i < NUMBER_OF_ARGS; ++i) {
      args[i] = atoi(argv[i + 1]);
    }
  }
  // Bcast args since MPI does not guarantee args to all PE:s
  MPI_Bcast(args, NUMBER_OF_ARGS, MPI_INT, ROOT, MPI_COMM_WORLD);

  // Get args unrolled
  N = args[0];
  iters = args[1];
  visualize = args[2];


  // Input arg for N x N per process
  const int block_sz = (N) / sqrt(comm_sz);

  // Calculate starting points:
  // Size of halo region is the block size padded with +2
  const int halo_sz = block_sz + 2;
  // Matrix size for halo matrix, for convenience
  const int halo_mat_sz = halo_sz * halo_sz;

  // Row dataype, excluding the corner, in contiguous memory
  MPI_Datatype row_type;
  MPI_Type_contiguous(block_sz, MPI_INT, &row_type);
  MPI_Type_commit(&row_type);

  // Column datatype containing the corner
  // which eliminates intercardinal communication
  // i.e diagonal exchange of corners.
  MPI_Datatype col_type;
  MPI_Type_vector(halo_sz, 1, halo_sz, MPI_INT, &col_type);
  MPI_Type_commit(&col_type);

  // Create dims for perfect squares
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  int dims[2] = {(int)sqrt(comm_sz), (int)sqrt(comm_sz)};
  MPI_Dims_create(comm_sz, 2, dims);
  int periods[2] = {true, true};
  int reorder = true;
  int nbrs[4];

  // Create a communicator with a cartesian topology.
  MPI_Comm cart_comm;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);
  MPI_Comm_rank(cart_comm, &cart_rank);

  // Find nearest neighbours
  MPI_Cart_shift(cart_comm, ROW_DIM, NEAREST_NBRS, &nbrs[0], &nbrs[1]);
  MPI_Cart_shift(cart_comm, COL_DIM, NEAREST_NBRS, &nbrs[2], &nbrs[3]);
  enum DIRECTIONS { WEST, EAST, SOUTH, NORTH };

  // Let MPI allocate row-major memory since it _might_ be faster in RMA.
  // The first half of the region contains 'next_state' followed
  // by 'this_state'
  MPI_Alloc_mem(2 * halo_sz * halo_sz * sizeof(int), MPI_INFO_NULL, &state);
  this_state = state + halo_mat_sz;
  next_state = state;

  // Initialize inner N x N matrix to random 0/1
  init_halo_mat(this_state, halo_sz);

  // Create window and find the subset of the window group
  // that contains the nearest neighbours ('win_nbrs')
  MPI_Group win_group, win_nbrs;
  MPI_Win window;
  MPI_Win_create(state, 2 * halo_sz * halo_sz * sizeof(int), sizeof(int),
                 MPI_INFO_NULL, cart_comm, &window);
  MPI_Win_get_group(window, &win_group);

  // Serial case where the root rank is its own neighbours
  if (comm_sz == 1) {
    MPI_Group_incl(win_group, comm_sz, nbrs, &win_nbrs);
  } else {
    MPI_Group_incl(win_group, NBRS, nbrs, &win_nbrs);
  }

  start = MPI_Wtime();
  for (int i = 0; i < iters; ++i) {
    // PSCW RMA epoch can be expressed in one block since all
    // PE:s belongs to every Put-operation. The offset allows
    // putting the boundary value alternating between
    // this_state and next_state.
    MPI_Win_post(win_nbrs, 0, window);
    MPI_Win_start(win_nbrs, 0, window);
    MPI_Put(&this_state[DOWN_ROW_POS], 1, row_type, nbrs[SOUTH],
            UP_HALO_POS + TIMESTEP_OFFSET, 1, row_type, window);
    MPI_Put(&this_state[UP_ROW_POS], 1, row_type, nbrs[NORTH],
            DOWN_HALO_POS + TIMESTEP_OFFSET, 1, row_type, window);
    MPI_Put(&this_state[RIGHT_COL_POS], 1, col_type, nbrs[EAST],
            LEFT_HALO_POS + TIMESTEP_OFFSET, 1, col_type, window);
    MPI_Put(&this_state[LEFT_COL_POS], 1, col_type, nbrs[WEST],
            RIGHT_HALO_POS + TIMESTEP_OFFSET, 1, col_type, window);
    MPI_Win_complete(window);
    MPI_Win_wait(window);

    // Sir Conway's Game of Life
    game_of_life(halo_sz, this_state, next_state);

    // Exchange pointers
    temp_state = next_state;
    next_state = this_state;
    this_state = temp_state;

    // Prints state from ROOT
    // Omitted in performance evaluation to avoid branching
    if (visualize == 1) {
      CLEAR_AND_PRINT_STATE
    }
  }

  stop = MPI_Wtime();

  time = stop - start;

  // Get the maximum runtime among the processes
  MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, ROOT, cart_comm);

  // Print reduced maximum
  if (cart_rank == ROOT) {
    printf("%d, %d, %d, %lf\n", iters, N, comm_sz, time_max);
  }
  
  // Free types, memory and finalize
  MPI_Type_free(&col_type);
  MPI_Type_free(&row_type);
  MPI_Win_free(&window);
  MPI_Free_mem(state);
  MPI_Finalize();
  return 0;
}
