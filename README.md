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

ctest

une série de tests démarre. Il est conseillé de faire repasser ces
tests après toute modification du code.

Pour lancer schnaps:

générer le maillage:

cd geo ; gmsh disque.geo -3 ; cd ..

./schnaps

Cet exemple consiste à résoudre l'équation de transport à
l'intérieur d'un cylindre. La visualisation  des résultats
utilise gmsh

gmsh dgvisu.msh

(puis tools -> options -> view[0] -> adapt visualization grid pour
afficher une image plus jolie...)

Adapter le fichier source "schnaps.c"
pour traiter des cas avec d'autres maillages. Des exemples de
maillages se trouvent dans le dossier geo.

Génération de la documentation (avec doxygen):

cd doc/
doxygen doxyschnaps


