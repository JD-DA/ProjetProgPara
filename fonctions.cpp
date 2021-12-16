//
// Created by sophie on 13/12/2021.
//

#include <mpi.h>
#include "Terrain.h"

void round_robin(mnt* mnt, float* terrain_local, int taille_bande, int nb_lignes, int nb_cols, int root) {
    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int nb_bandes = (nb_lignes / taille_bande) / (nprocs);
    int reste = (nb_lignes / taille_bande) % nprocs;

    int ptr_terrain;
    int n_envoi = (taille_bande+2)*nb_cols;
    int n_envoi_reduit = (taille_bande+1)*nb_cols;
    if (pid == root) {
        float *terrain = mnt->get_data();
        // Première bande avec une bordure en moins pour le processus de pid 0
        if (root != 0) {
            ptr_terrain = 0;
            MPI_Ssend(terrain + ptr_terrain, n_envoi_reduit, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        } else {
            for (int l = 0; l < n_envoi_reduit; l++) {
                terrain_local[nb_cols + l] = terrain[l];
            }
        }
        // Première bande pour les autres processus
        for (int dest = 1; dest < nprocs; dest++) {
            ptr_terrain=(dest*taille_bande-1)*nb_cols;
            if (dest != root) {
                MPI_Ssend(terrain +ptr_terrain, n_envoi, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
            } else {
                for (int l = 0; l < n_envoi; l++)
                    terrain_local[l] = terrain[ptr_terrain + l];
            }
        }
        // Distribution des bandes (ni la première ni la dernière) le cas général
        for (int i = 1; i < nb_bandes-1; i++) {
            for (int dest = 0; dest < nprocs; dest++) {
                ptr_terrain = ((i*nprocs + dest)*taille_bande-1)*nb_cols;
                if (dest != root) {
                    MPI_Ssend(terrain + ptr_terrain, n_envoi, MPI_FLOAT, dest, i, MPI_COMM_WORLD);
                }
                else {
                    for (int l=0; l<(taille_bande+2)*nb_cols; l++)
                        terrain_local[i*(taille_bande+2)*nb_cols+l] = terrain[ptr_terrain+l];
                }
            }
        }
        // Dernière bande avec une bordure en moins pour le processus de pid nprocs-1
        if (root!=nprocs-1) {
            ptr_terrain=(((nb_bandes-1)*nprocs + nprocs-1)*taille_bande-1)*nb_cols;
            MPI_Ssend(terrain+ptr_terrain,n_envoi_reduit, MPI_FLOAT, nprocs-1, nb_bandes-1, MPI_COMM_WORLD);
        }
        else {
            for (int l = 0; l < (taille_bande + 1) * nb_cols; l++) {
                terrain_local[(nb_bandes-1)*(taille_bande+2)*nb_cols+l] = terrain[ptr_terrain+l];
            }
        }
        // Dernière bande pour les autres processus
        for (int dest=0; dest<nprocs-1; dest++) {
            ptr_terrain=(((nb_bandes-1)*nprocs + dest)*taille_bande-1)*nb_cols;
            if (dest!=root){
                MPI_Ssend(terrain+ptr_terrain,n_envoi, MPI_FLOAT, dest, nb_bandes-1, MPI_COMM_WORLD);
            }
            else {
                for (int l = 0; l < (taille_bande + 2) * nb_cols; l++) {
                    terrain_local[(nb_bandes-1)*(taille_bande+2)*nb_cols+l] = terrain[ptr_terrain+l];
                }
            }
        }
    }
    else { // Pour la partie réception
        if (pid == 0) // la première bande sur le processus pid 0
            MPI_Recv(terrain_local + nb_cols, n_envoi_reduit, MPI_FLOAT, root, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        else { // la première bande sur les autres processus
            MPI_Recv(terrain_local, n_envoi, MPI_FLOAT, root, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }
        for (int i = 1; i < nb_bandes - 1; i++) { // Ici le cas général
            MPI_Recv(terrain_local + i * (taille_bande + 2) * nb_cols, (taille_bande + 2) * nb_cols, MPI_FLOAT, root, i,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (pid == nprocs-1) { // la dernière bande sur le processus pid nproc -1
            MPI_Recv(terrain_local + (nb_bandes - 1) * (taille_bande + 2) * nb_cols, (taille_bande + 1) * nb_cols,
                     MPI_FLOAT, root, nb_bandes-1, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }
        else { // la dernière bande pour les autres processus.
            MPI_Recv(terrain_local + (nb_bandes - 1) * (taille_bande + 2) * nb_cols, (taille_bande + 2) * nb_cols,
                     MPI_FLOAT, root, nb_bandes-1, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }

    }
}

void calcul_direction(float *terrain_local, int *dir, int nb_bandes, int taille_bande, int nb_cols, float no_value){
    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int nb_bande_local = nb_bandes/nprocs;
    for (int i = 0; i < nb_bande_local; ++i) {
#pragma omp parallel for
        for (int y = 0; y < taille_bande; ++y) {
            for (int x = 0; x < nb_cols; ++x) {
                if(x==0){
                    dir[coords_to_indice(x, y, nb_bande_local, taille_bande, nb_cols)] = chercher_min_bord(x, y);
                }else if(x==nb_cols-1){
                    dir[coords_to_indice(x, y, nb_bande_local, taille_bande, nb_cols)] = chercher_min_bord(x, y);
                }else if(y==0 and pid==0 and nb_bande_local==0){
                    dir[coords_to_indice(x, y, nb_bande_local, taille_bande, nb_cols)] = chercher_min_bord(x, y);
                }else if (y==taille_bande-1 and pid==nprocs-1 and nb_bande_local==nb_bandes-1){
                    dir[coords_to_indice(x, y, nb_bande_local, taille_bande, nb_cols)] = chercher_min_bord(x, y);
                }else{
                    dir[coords_to_indice(x, y, nb_bande_local, taille_bande, nb_cols)] = chercher_min(x, y);
                }
            }
        }
    }
}

int coords_to_indice(int x,int y, int nBande,int taille_bande,int nb_cols){
    return x+nb_cols*(1+y+nBande*(taille_bande+2));
}

int chercher_min(int x, int y){
    return 1;
}

int chercher_min_bord(int x, int y){
    return 1;
}


