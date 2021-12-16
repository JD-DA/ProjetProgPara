#include <iostream>
#include <mpi.h>
#include "Terrain.h"

int main(int argc, char **argv) {
    int pid, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   
    char* filename = argv[1];
    int root = atoi(argv[2]);
    mnt* donnees;

    if (pid==root) {
        donnees = new mnt(filename);
        cout << *donnees << endl;
    }

    if (pid==root)
      delete donnees;
    MPI_Finalize();
    return 0;
}
