//
// Created by o2184043 on 16/12/2021.
//

#include <iostream>
#include "mnt.h"
#include "Terrain.h"

void round_robin(mnt* mnt, float* terrain_local, int taille_bande, int nb_lignes, int nb_cols, int root) {

    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nbprocs);

    int nb_bandes = nb_lignes / taille_bande;
    int taille_loc = (taille_bande+2)*nb_cols*nb_bandes/nprocs;

    if (nb_lignes % taille_bande != 0)
        return;

    if (nb_bandes % nprocs != 0)
        return;

    float *tab = new float [taille_loc];

    int debut_rcv = 1*nb_cols;
    int debut_envoi;

    //Envois des bandes aux processeurs
    for (int i=0; i<nb_bandes/nprocs; i++) {
        debut_envoi = taille_bande*nprocs*i*nb_cols;
        MPI_Scatter(mnt->get_data()+debut_envoi, taille_bande*nb_cols, MPI_FLOAT,
                tab+debut_rcv, taille_bande*nb_cols, MPI_FLOAT,
                root, MPI_COMM_WORLD);
        debut_rcv += (taille_bande+2)*nb_cols;
    }

    //Echange des lignes voisines
    int voisin_haut = (pid-1)%nprocs;
    int voisin_bas = (pid+1)%nprocs;

    int debut_rcv_haut = 0;
    int debut_envoi_haut = 1*nb_colsb;

    int debut_rcv_bas = 0;
    int debut_envoi_bas = 1*nb_colsb;

    for (int i=0; i<nb_bandes/nprocs; i++) {
        if (pid % 2 == 0) {

            MPI_Send(tab+nb_cols, nb_cols, MPI_FLOAT, voisin_haut, 1, MPI_COMM_WORLD);
            MPI_Recv(tab, nb_cols, MPI_FLOAT, voisin_haut, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(tab+nb_cols, nb_cols, MPI_FLOAT, voisin_bas, 1, MPI_COMM_WORLD);
            MPI_Recv(tab, nb_cols, MPI_FLOAT, voisin_bas, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else {

            MPI_Send(tab+nb_cols, nb_cols, MPI_FLOAT, voisin_haut, 1, MPI_COMM_WORLD)
        }
    }
}