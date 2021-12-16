#Partie 3 : Directions
###Travail réalisé :
Nous avons réalisé une fonction qui va itérer sur chaque bande locale puis l'ensemble des cases de cette bande afin de chercher le minimum autour de celle-ci et donc la direction que l'on doit lui attribuer.
De plus pour permettre une certaine parallélisation nous avons ajouté une commande ```#prama omp parallel for``` au niveau du ```for``` qui itére sur nos lignes de bande.
Afin d'utiliser nos bandes et les ghosts disponible nous avons implémenté une fonction ```coords_to_indice``` qui pour un couple x,y donné retourne l'emplacement dans ```terrain_local``` de cette donnée.
Cette fonction prend en compte les plage de données ghosts et renvoi aussi leurs positions.

###Justification du ```#prama omp parallel for```
ON se trouve sur une machine hybride, via MPI on parallélise les calculs en répartissant les données sur les différents noeuds de la machine mais afin d'optimiser un peu plus notre programme on peut aussi dans les calculs sur chaque noeud de rajouter une directive omp qui va donc répartir la charge du for sur les différents coeurs du processeur du noeud de la machine. 
Ainsi ici pour chaque bande on va répartir les lignes calculées.
De plus on peut utiliser un ```schedule(dynamic)```


#Partie 4 : Accumulation
###Travail réalisé :
Pour cette partie, nous avons écrit une fonction qui permet de communiquer les résultats du calcul des directions aux voisins.
Nous avons également écrit la fonction qui permet de calculer les accumulations. Cette dernière utilise une variante de la fonction ci-dessus pour communiquer les ghost aux bandes voisines entre chaque itération.
Une itération de la fonction consiste à parcourir tous les points pas encore marqués et vérifier s'il est possible de calculer leur flot (c'est à dire que les voisins qui se déverse sur le point sont déjà marqués).
