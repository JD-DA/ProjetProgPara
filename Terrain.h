//
// Created by sophie on 10/12/2021.
//
#include <iostream>
using namespace std;
class mnt {
private:
    float *terrain;
    int n_lignes;
    int n_cols;
    float no_value;

public:
    mnt(char *nom);
    ~mnt();

    int get_nb_lignes();

    int get_nb_cols();

    float* get_data();

    float get_no_value();

  friend ostream& operator<<(ostream& ost, const mnt& data);
};
