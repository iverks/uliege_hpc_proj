# HPC Project

## Important commands 

`data2gmsh -cutx 0.5 ./out_cutx_p.dat`

Remember to set OMP_NUM_THREADS on nic5 jobs using `export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK`

## Analysis with scalesca

Score-p to add probe points `scorep mpicc -O3 -fopenmp -o ...`

Run with scalasca run `scalasca -analyze srun ../../fdtd param3d.txt`
Examine with scalasca `scalasca -examine -f ./filter_file.filter  ./mappe_navn`, output is a file named `something.score`

Run with scalasca run disable profiling and enable tracing:

```bash
export SCOREP_TOTAL_MEMORY=256MB

scalesca -analyze -q -t -f filter_file.filter srun ./myfile
```

scp the `something.cubex` file and open in cubex

## Analysis with LIKWID

In the hardware there are counters, we can access these with likwid. (perf as well but it is more complex, only for nerds)



