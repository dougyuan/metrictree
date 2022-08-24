# Parallelized Metric Tree Construction #
Final project for CS 140: Parallel Scientific Computing

# How to Use #
To compile all three versions of the program, use:

> $ make all

Benchmarks can be performed by using:
> $ ./Version number_of_points dimensions

where
- Version is any of the following: <i>serial, pthreads, openmp</i>.
- Point and dimension values must be positive integers.

Example:
> $ ./openmp 1000000 50
