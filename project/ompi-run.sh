mpicc tsp_omp.c -o tsp_omp -fopenmp -lm -ldl
mpirun -n 5 ./tsp_omp

