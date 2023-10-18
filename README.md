# HPC Project

## Important commands 

`data2gmsh -cutx 0.5 ./out_cutx_p.dat`

Remember to set OMP_NUM_THREADS on nic5 jobs using `export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK`

## Task 1

This was pretty straightforward, just made a function called `double trilinear_interpolate(data_t dat, double x, double y, double z)`
that interpolates over the grid of `data_t` and gets the proper approximation. Very manually coded, nothing super clever.

## Task 2
 
$$
\frac{c \Delta_t }{\delta}
$$

Varying $\Delta_t$. $\delta = 1e-2$, $c = const$

| Delta_t | t_tot | comment           |
| ------- | ----- | ----------------- |
| 5e-7    | 0.001 | stable            |
| 1e-6    | 0.002 | relatively stable |
| 1e-5    | 0.02  | unstable          |
| 4e-6    | 0.008 | stable            |

## Task 3

Per loop: 

- Update pressure
  - read 9 values
    - each indexing uses 3 multiplications and 2 additions
  - compute 3 + 1 multiplications, 3 + 1 subtraction, 2 additions, 3 comparisons


## Task 4

This was solved by changing the index3D macro from

```c
#define INDEX3D(grid, m, n, p) \
    ((size_t)grid.numnodesy * grid.numnodesx * (p) + grid.numnodesx * (n) + (m))
```

to 

```c
#define INDEX3D(grid, m, n, p) \
    ((size_t)(p) + grid.numnodesz * (n) + grid.numnodesz * grid.numnodesy * (m))
```

which changes the matrix from having values with the same p (z) and n (y), but different m (x) after each 
other to having values with the same m (x) and n (y), but different p (z), which is important because the
innermost loop in the code is the looping over p, so we want values with p +/- 1 to 
be in the same cache line.