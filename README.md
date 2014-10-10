SCHNAPS
=======

Solveur pour les lois de Conservation Hyperboliques Non-linéaires Appliqué aux PlasmaS

Mode d'emploi

téléchargement (nécessite git)

git clone https://github.com/phelluy/SCHNAPS.git

se placer dans le dossier SCHNAPS

puis:

cmake .

make

./schnaps

une série de tests démarre. Il est conseillé de faire repasser ces
tests après toute modification du code.

Le dernier test consiste à résoudre l'équation de transport à
l'intérieur d'un cylindre. La visualisation  des résultats
utilise gmsh

gmsh dgvisu.msh

(puis tools -> options -> view[0] -> adapt visualization grid pour
afficher une image plus jolie...)

Adapter la fonction "TestFieldRK2" dans le fichier source "test.c"
pour traiter des cas avec d'autres maillages. Des exemples de
maillages se trouvent dans le dossier geo.
Pour générer un maillage à partir du fichier disque.geo par exemple,
faire:

gsmh disque.geo -3

Cette commande génère un fichier disque.msh lisible par SCHNAPS.


