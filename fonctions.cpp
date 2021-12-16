//
// Created by sophie on 13/12/2021.
//

#include <mpi.h>
#include "Terrain.h"
#include "fonctions.h"
#include <queue>

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
        if(pid==0){
            cout<<i<<endl;
        }
#pragma omp parallel for schedule(runtime)
        for (int y = 0; y < taille_bande; ++y) {
            for (int x = 0; x < nb_cols; ++x) {

                if(x==0){
                    dir[coords_to_indice(x, y, i, taille_bande, nb_cols)] = chercher_min_bord(x, y);
                }else if(x==nb_cols-1){
                    dir[coords_to_indice(x, y, i, taille_bande, nb_cols)] = chercher_min_bord(x, y);
                }else if(y==0 and pid==0 and i==0){
                    dir[coords_to_indice(x, y, i, taille_bande, nb_cols)] = chercher_min_bord(x, y);
                }else if (y==taille_bande-1 and pid==nprocs-1 and i==nb_bandes-1){
                    dir[coords_to_indice(x, y, i, taille_bande, nb_cols)] = chercher_min_bord(x, y);
                }else{
                    dir[coords_to_indice(x, y, i, taille_bande, nb_cols)] = chercher_min(x, y,terrain_local,i,taille_bande,nb_cols,no_value);
                }
            }
        }
    }
}

int coords_to_indice(int x,int y, int nBande,int taille_bande,int nb_cols){
    return x+nb_cols*(1+y+nBande*(taille_bande+2));
}

/*
 * Va chercher le min autour du point de coordonnées x,y en tenant compte de la possiblité d'un nodata
 */
int chercher_min(int x, int y,float *terrain_local, int nBande, int taille_bande, int nb_cols,float nodata){
    float tab[9];
    if(terrain_local[coords_to_indice(x,y,nBande,taille_bande,nb_cols)]==nodata)
        return 0;
    tab[0]=terrain_local[coords_to_indice(x,y,nBande,taille_bande,nb_cols)];
    tab[1]=terrain_local[coords_to_indice(x,y-1,nBande,taille_bande,nb_cols)];
    tab[2]=terrain_local[coords_to_indice(x+1,y-1,nBande,taille_bande,nb_cols)];
    tab[3]=terrain_local[coords_to_indice(x+1,y,nBande,taille_bande,nb_cols)];
    tab[4]=terrain_local[coords_to_indice(x+1,y+1,nBande,taille_bande,nb_cols)];
    tab[5]=terrain_local[coords_to_indice(x,y+1,nBande,taille_bande,nb_cols)];
    tab[6]=terrain_local[coords_to_indice(x-1,y+1,nBande,taille_bande,nb_cols)];
    tab[7]=terrain_local[coords_to_indice(x-1,y,nBande,taille_bande,nb_cols)];
    tab[8]=terrain_local[coords_to_indice(x-1,y-1,nBande,taille_bande,nb_cols)];
    float mini=tab[0];
    int imin=0;
    for (int i = 1; i < 9; ++i) {
        if(tab[i]<mini and tab[i]!=nodata){
            mini=tab[i];
            imin=i;
        }
    }
    return imin;

}

/**
 * Fonction spéciale qui va chercher le min dans le cas où l'on se trouve sur un bord ou un angle. TO DO
 * @param x
 * @param y
 * @return
 */
int chercher_min_bord(int x, int y){
    return 1;
}


void calcul_accumulation(float *terrain_local, int *dir, float *acc, int nb_bandes, int taille_bande, int nb_cols, float no_value) {
    int pid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int nb_bandes_local = nb_bandes/nprocs;

    //Initialisation
    int nb_non_marques = 0;
    queue<int*> points_non_marques;
    for (int i=0; i<nb_bandes_local; ++i) {
        for (int y=0; y<taille_bande; ++y) {
            for (int x=0; x<nb_cols; ++x) {
                int j = coords_to_indice(x,y,i,taille_bande,nb_cols);
                if (dir[j] == 0 or terrain_local[j] == no_value)
                    acc[j] = 1.;
                else {
                    acc[j] = 0.;
                    nb_non_marques++;
                    int triplet[3] = {x,y,i};
                    points_non_marques.push(triplet); //Un point de la file est définit par des coordonées dans une bande et le numéro de cette bande
                }
            }
        }
    }

    while (nb_non_marques > 0) { //On continue de boucler tant qu'il reste des points à calculer
        queue<int*> new_points_non_marques;
        while(!points_non_marques.empty()) { //On parcourt la file des points non-marqués pour tester s'il sont calculables
            int *coords = points_non_marques.front();
            points_non_marques.pop();
            int x = coords[0];
            int y = coords[1];
            int num_bande = coords[2];
            float *res_voisins = verifier_voisins(x,y,acc,dir,num_bande,taille_bande,nb_cols);
            float somme = 1.;
            bool calculer_point = true;
            for (int j=1; j<9; ++j) { //Quand on parcourt les voisins, si une des valeurs est égale à 0, cela veut dire qu'il est dirigé vers le point mais n'est pas encore marqué
                float res_j = res_voisins[j];
                if (res_j == 0.)
                    calculer_point = false;
                else if (res_j > 0.)
                    somme += res_j;
            }
            if (calculer_point) {
                acc[coords_to_indice(x, y, num_bande, taille_bande, nb_cols)] = somme;
                nb_non_marques--;
            } else {
                new_points_non_marques.push(coords); //Si l'on n'a pas calculé le point, on le remet dans la nouvelle file
            }
        }
        points_non_marques = new_points_non_marques;

        //Comunications des ghost entre bandes voisines
    }
}

/**
 * @param x
 * @param y
 * @param acc
 * @param dir
 * @param nBande
 * @param taille_bande
 * @param nb_cols
 * @return un tableau de float où, pour le voisin i, tab[i] est -1 s'il ne se déverse pas sur le point, et sinon tab[i] est la valeur de son flot
 */
float* verifier_voisins(int x, int y, float *acc, int *dir, int nBande, int taille_bande, int nb_cols) {
    float *resultats_voisins = new float [9];

    int j = coords_to_indice(x-1,y-1,nBande,taille_bande,nb_cols); //point nord-ouest
    if (dir[j] == 4)
        resultats_voisins[8] = acc[j];
    else
        resultats_voisins[8] = -1.; //Valeur qui sera ignorée dans le parcours des voisins

    j = coords_to_indice(x,y-1,nBande,taille_bande,nb_cols); //point nord
    if (dir[j] == 5)
        resultats_voisins[1] = acc[j];
    else
        resultats_voisins[1] = -1.;

    j = coords_to_indice(x+1,y-1,nBande,taille_bande,nb_cols); //point nord-est
    if (dir[j] == 6)
        resultats_voisins[2] = acc[j];
    else
        resultats_voisins[2] = -1.;

    j = coords_to_indice(x-1,y,nBande,taille_bande,nb_cols); //point ouest
    if (dir[j] == 3)
        resultats_voisins[7] = acc[j];
    else
        resultats_voisins[7] = -1.;

    j = coords_to_indice(x+1,y,nBande,taille_bande,nb_cols); //point est
    if (dir[j] == 7)
        resultats_voisins[3] = acc[j];
    else
        resultats_voisins[3] = -1.;

    j = coords_to_indice(x-1,y+1,nBande,taille_bande,nb_cols); //point sud-ouest
    if (dir[j] == 2)
        resultats_voisins[6] = acc[j];
    else
        resultats_voisins[6] = -1.;

    j = coords_to_indice(x,y+1,nBande,taille_bande,nb_cols); //point sud
    if (dir[j] == 1)
        resultats_voisins[5] = acc[j];
    else
        resultats_voisins[5] = -1.;

    j = coords_to_indice(x+1,y+1,nBande,taille_bande,nb_cols); //point sud-est
    if (dir[j] == 8)
        resultats_voisins[4] = acc[j];
    else
        resultats_voisins[4] = -1.;

    return resultats_voisins;
}