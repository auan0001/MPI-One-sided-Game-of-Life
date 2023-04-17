# Game of Life One-sided MPI

## About

This is an excerpt from my BSc thesis in the Engineering Physics programme at the dept. of Computational Science, Uppsala universitet. From this repository, the One-sided implementation as well as a two-sided `MPI_Sendrecv` (for comparison) can be compiled to run on a system with an installation of `OpenMPI`. The data passing between the nearest neighbors is done through halo exchange by using derived datatypes of boundary rows and columns

![Halo exchange](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_halo.png)

## Compiling and running
Clone the repository and run `Make` in the root directory.

Running the program

```
mpirun -np <NPROCS> <ONESIDED_GOL> <SIZE> <ITERS> <VISUALIZE (1 or 0)>
```

Where the visualization shows one arbitrary rank for debugging purposes. If the grid dimensions are incorrect


```
Either <NPROCS> not a perfect square or \n<PROBLEM SIZE> not divisible by sqrt<NPROCS>
```

will be written to the console.

## Results
This is results from running the code on the [UPPMAX Rackham](https://www.uppmax.uu.se/resources/systems/the-rackham-cluster/) cluster using the modules `GCC 9.2.0` and `OpenMPI 4.0.2`.

### Runtime
![Runtime](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_run.png)

### Strong scaling
![Strong scaling](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_str.png)

### Weak scaling
![Weak scaling](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_weak.png)

### Efficiency
![Efficiency](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_eff.png)

## Acknowledgements
The computations was enabled by resources in project snic2021-22-633 provided by the Swedish National Infrastructure for Computing (SNIC) at UPPMAX, partially funded by the Swedish Research Council through grant agreement no. 2018-05973.
