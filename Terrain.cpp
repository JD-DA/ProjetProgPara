//
// Created by sophie on 10/12/2021.
//
#include <iostream>
#include "Terrain.h"

mnt::mnt(char *nom) {
    FILE *f;
    f = fopen(nom, "r");
    if (f != NULL) {
        int tmp;
        fscanf(f, "%d", &tmp);
        this->n_lignes = tmp;
        fscanf(f, "%d", &tmp);
        this->n_cols = tmp;

        float novalue;

        fscanf(f, "%d", &tmp);
        fscanf(f, "%d", &tmp);
        fscanf(f, "%d", &tmp);
        fscanf(f, "%f", &novalue);
        this->no_value = novalue;

        this->terrain = new float[this->n_cols*this->n_lignes];

        for (int i = 0; i < this->n_lignes; i++)
            for (int j = 0; j < this->n_cols; j++)
                fscanf(f, "%f", &(this->terrain[i * this->n_cols + j]));
    }
}

mnt::~mnt()
{
    delete[] terrain;
}
int mnt::get_nb_lignes(){
    return n_lignes;
};

    int mnt::get_nb_cols(){
        return n_cols;
    }
float* mnt::get_data(){
        return terrain;
    }
float mnt::get_no_value(){
        return no_value;
    }
ostream& operator<<(ostream& ost, const mnt& data)
{
    for (int i=0; i<data.n_lignes; i++) {
        for (int j=0; j<data.n_cols; j++)
            ost << data.terrain[i*data.n_cols+j] << " ";
        ost << endl;
    }
    return ost;
}