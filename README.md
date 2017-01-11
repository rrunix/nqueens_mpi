# nqueens_mpi

Nqueens MPI version (using the sequential solution proposed by Jeff Somers (http://www.jsomers.com/nqueen_demo/nqueens.html)).

Usage:
```
  mpicc *.c -o nq -lm
  
  ./nq <width of board> <server depth>
```
Where server depth is the number column number from which the master delegates to the slaves, ie, the master only do the first
server depth rows and then left the rest for the slaves.
