//
// Created by o2184043 on 16/12/2021.
//

#ifndef PROJETPROGPARA_MNT_H
#define PROJETPROGPARA_MNT_H

void round_robin(mnt* mnt, float* terrain_local, int taille_bande, int nb_lignes, int nb_cols, int root);

void calcul_direction(float *terrain_local, int *dir, int nb_bandes, int taille_bande, int nb_cols, float no_value);

#endif //PROJETPROGPARA_MNT_H
