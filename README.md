# Game of Life One-sided MPI
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
This is results from running the code on the [UPPMAX Rackham](https://www.uppmax.uu.se/resources/systems/the-rackham-cluster/) cluster.

### Runtime
![Runtime](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_run.png)

### Strong scaling
![Strong scaling](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_str.png)

### Weak scaling
![Weak scaling](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_weak.png)

### Efficiency
![Efficiency](https://github.com/auan0001/MPI-One-sided-Game-of-Life/blob/main/images/gol_gh_eff.png)
