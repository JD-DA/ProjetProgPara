//
// Created by sophie on 13/12/2021.
//

#ifndef MNT_DIRECTION_ET_ACCUMULATION_FONCTIONS_H
#define MNT_DIRECTION_ET_ACCUMULATION_FONCTIONS_H

void round_robin(mnt* mnt, float* terrain_local, int taille_bande, int nb_lignes, int nb_cols, int root);
void calcul_direction(float *terrain_local, int *dir, int nb_bandes, int taille_bande, int nb_cols, float no_value);
void calcul_accumulation(float *terrain_local, float *dir, float *acc, int nb_bandes, int taille_bande, int nb_cols, float no_value);

//fonctions auxiliaires
int coords_to_indice(int x,int y, int nBande,int taille_bande,int nb_cols);
int chercher_min_bord(int x, int y);
int chercher_min(int x, int y,float *terrain_local, int nBande, int taille_bande, int nb_cols,float nodata);

float* verifier_voisins(int x, int y, float *acc, int *dir, int nBande, int taille_bande, int nb_cols);
#endif //MNT_DIRECTION_ET_ACCUMULATION_FONCTIONS_H
