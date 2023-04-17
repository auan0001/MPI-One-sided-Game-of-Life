#include "utils.h"

int main(int argc, char *argv[]) {
  // Communicator size and ranks
  int comm_sz, rank, cart_rank;
  // Input args
  int args[3], N, iters, visualize;
  // Memory pointers
  int *state, *this_state, *temp_state, *next_state;
  double start, stop, time, time_max;

  MPI_Init(&argc, &argv);

  // Size of the default communicator
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Random seed dependent on rank
  srand(rank*1337);

  // Check input args
  if (rank == ROOT) {
    for (int i = 0; i < 3; ++i) {
      args[i] = atoi(argv[i+1]);
    }
    if (argc < 4) {
      printf("--------------------------------------------------------------------------\n");
      printf("Usage: mpirun -np <NPROCS> <ONESIDED_GOL> <SIZE> <ITER> <VISUALIZE>\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (comm_sz%(int)sqrt(comm_sz) || atoi(argv[1])%((int)sqrt(comm_sz))) {
      printf("--------------------------------------------------------------------------\n");
      printf("Either <NPROCS> not a perfect square or \n<PROBLEM SIZE> not divisible by sqrt<NPROCS>\n");
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
  // Bcast argv since MPI does not guarantee args to all PE:s
  MPI_Bcast(args, 3, MPI_INT, ROOT, MPI_COMM_WORLD);

  N = args[0];
  iters = args[1];
  visualize = args[2];

  // Input arg for N x N per process
  const int block_sz = (N)/sqrt(comm_sz); 

  // Calculate starts
  const int halo_sz = block_sz+2;
  const int halo_mat = halo_sz*halo_sz;

  // Row dataype in contiguous memory
  MPI_Datatype row_type;
  MPI_Type_contiguous(block_sz, MPI_INT, &row_type);
  MPI_Type_commit(&row_type);

  // Column datatype containing the corner
  // which eliminates intercardinal communication
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
  enum DIRECTIONS {WEST, EAST, SOUTH, NORTH};

  // Let MPI allocate row-major memory since it _might_ be faster in RMA.
  // The first half of the region contains 'next_state' followed 
  // by 'this_state'
  MPI_Alloc_mem(2*halo_sz*halo_sz*sizeof(int), MPI_INFO_NULL, &state);
  this_state = state+halo_mat;
  next_state = state;

  // Initialize inner N x N matrix to random 0/1
  init_halo_mat(this_state, halo_sz);

  start = MPI_Wtime();
  for (int i = 0; i < iters; ++i) {

    // Send and receive within the cartesian neighbours
    MPI_Sendrecv(&this_state[UP_ROW_POS], 1, row_type, nbrs[SOUTH],
        TAG, &state[DOWN_HALO_POS+TIMESTEP_OFFSET], 1, row_type, nbrs[NORTH],
        TAG, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&this_state[DOWN_ROW_POS], 1, row_type, nbrs[NORTH],
        TAG, &state[UP_HALO_POS+TIMESTEP_OFFSET], 1, row_type, nbrs[SOUTH],
        TAG, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&this_state[LEFT_COL_POS], 1, col_type, nbrs[EAST],
        TAG, &state[RIGHT_HALO_POS+TIMESTEP_OFFSET], 1, col_type, nbrs[WEST],
        TAG, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&this_state[RIGHT_COL_POS], 1, col_type, nbrs[WEST],
        TAG, &state[LEFT_HALO_POS+TIMESTEP_OFFSET], 1, col_type, nbrs[EAST],
        TAG, cart_comm, MPI_STATUS_IGNORE);

    // Sir Conway's Game of Life
    game_of_life(halo_sz, this_state, next_state);

    // Exchange pointers
    temp_state = next_state;
    next_state = this_state;
    this_state = temp_state;

    if (visualize == 1) {CLEAR_AND_PRINT_STATE}
  }

  stop = MPI_Wtime();

  time = stop-start;

  // Get the maximum runtime from one of the processes
  MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, ROOT, cart_comm);

  if (cart_rank == ROOT){
    printf("%d, %d, %d, %lf\n", iters, N, comm_sz, time_max);
  }

  MPI_Type_free(&col_type);
  MPI_Type_free(&row_type);
  MPI_Free_mem(state);
  MPI_Finalize();
  return 0;
}
