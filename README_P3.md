#Partie 3 : Directions
###Travail réalisé :
Nous avons réalisé une fonction qui va itérer sur chaque bande locale puis l'ensemble des cases de cette bande afin de chercher le minimum autour de celle-ci et donc la direction que l'on doit lui attribuer.
De plus pour permettre une certaine parallélisation nous avons ajouté une commande ```#prama omp parallel for``` au niveau du ```for``` qui itére sur nos lignes de bande.
Afin d'utiliser nos bandes et les ghosts disponible nous avons implémenté une fonction ```coords_to_indice``` qui pour un couple x,y donné retourne l'emplacement dans ```terrain_local``` de cette donnée.
Cette fonction prend en compte les plage de données ghosts et renvoi aussi leurs positions.