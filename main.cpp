#include <iostream>
#include <mpi.h>
#include <omp.h>
#include "Terrain.h"
#include "fonctions.h"

int main(int argc, char **argv) {
    int pid, nprocs;
    int provided;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    char* filename = argv[1];
    int root = atoi(argv[2]);
    int taille_bande =atoi(argv[3]);

    mnt* donnees;

    if (pid==root) {
        donnees = new mnt(filename);
       // cout << *donnees << endl;
    }
    int taille[2];
    if (pid==root) {
        taille[0] = donnees->get_nb_lignes();
        taille[1] = donnees->get_nb_cols();
    }
    MPI_Bcast(taille,2,MPI_INT,root,MPI_COMM_WORLD);
    float no_value;
    if (pid==root)
        no_value=donnees->get_no_value();

    MPI_Bcast(&no_value,1,MPI_FLOAT,root,MPI_COMM_WORLD);

    int nb_lignes = taille[0];
    int nb_cols = taille[1];
    int nb_bandes = (taille[0]/taille_bande)/(nprocs);

    cout << "je suis : " << pid << " taille : " << taille[0] << " " << taille[1] << " nb bandes : " << nb_bandes << endl;
    float* terrain_local = new float[(taille_bande + 2) * nb_bandes * nb_cols];
    std::fill_n(terrain_local, (taille_bande + 2) * nb_bandes * nb_cols,no_value);

    round_robin(donnees,terrain_local,taille_bande,nb_lignes,nb_cols,root);

    int* terrain_dir = new int[(taille_bande + 2) * nb_bandes * nb_cols];
    calcul_direction(terrain_local, terrain_dir, nb_bandes,taille_bande,nb_cols,no_value);
    //rassembler(terrain_local, terrain_dir,nb_bandes,taille_bande,nb_cols,no_value);
    /*cout<<"Tab directions local : "<<endl;
    if(pid==root){
        for (int i = 0; i < (taille_bande + 2) * nb_bandes * nb_cols; ++i) {
            if(i%nb_cols==0 and i!=0)
                cout<<endl;
            cout << terrain_dir[i] <<" ";
        }
        cout<<endl;
    }*/


   // cout << "je suis : " << pid << " taille : " << taille[0] << " " << taille[1] << " nb bandes : " << nb_bandes << endl;

    if (pid==1) {
       cout << "je suis " << pid << " terrain_local:";
       /*for (int i = 0; i < nb_bandes * (taille_bande + 2); i++)
           for (int j = 0; j < nb_cols; j++)
               cout << terrain_local[i * nb_cols + j] << " ";
       cout << endl;*/
   }
    MPI_Finalize();
    return 0;
}
