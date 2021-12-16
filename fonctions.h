//
// Created by sophie on 13/12/2021.
//

#ifndef MNT_DIRECTION_ET_ACCUMULATION_FONCTIONS_H
#define MNT_DIRECTION_ET_ACCUMULATION_FONCTIONS_H

void round_robin(mnt* mnt, float* terrain_local, int taille_bande, int nb_lignes, int nb_cols, int root);
void calcul_direction(float *terrain_local, int *dir, int nb_bandes, int taille_bande, int nb_cols, float no_value);
int coords_to_indice(int x,int y, int nBande,int taille_bande,int nb_cols);
#endif //MNT_DIRECTION_ET_ACCUMULATION_FONCTIONS_H
