//! \file schnaps.h Main header of the SCHNAPS library. Has to be included
//! by a program using the library.
#ifndef _SCHNAPS_H
#define _SCHNAPS_H



/*! \mainpage  Solveur pour les lois de Conservation Hyperboliques Non-linéaires Appliqué aux PlasmaS


\section quoi Quesaco ?
C'est un solveur pour résoudre des systèmes hyperboliques de lois de conservations. Il est basé sur la méthode Galerkin Discontinu. SCHNAPS peut exploiter les accélérateurs type GPU via OpenCL.

\section lic Licence

SCHNAPS est sous licence CeCILL:

http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.html

(very similar to GPL)


\section donload Téléchargement

nécessite git

Accès développeur:
git clone git+ssh://<gforge_account_name>\@scm.gforge.inria.fr//gitroot/schnaps/schnaps.git

Accès lecture seule:
git clone https://gforge.inria.fr/git/schnaps/schnaps.git

\section compil Compilation, tests

se placer dans le dossier SCHNAPS

puis:

cmake . (ne pas oublier le point ".")

make

ctest

une série de tests démarre. Il est conseillé de faire repasser ces
tests après toute modification du code.

Pour lancer schnaps:

générer le maillage:

gmsh disque.geo -3

./schnaps

Cet exemple consiste à résoudre l'équation de transport à
l'intérieur d'un cylindre. La visualisation  des résultats
utilise gmsh

gmsh dgvisu.msh

(puis tools -> options -> view[0] -> adapt visualization grid pour
afficher une image plus jolie...)

Adapter le fichier source "schnaps.c" pour traiter des cas avec
d'autres maillages. Des exemples de maillages se trouvent dans le
dossier geo.


\section doc Documentation

Génération de la documentation (avec doxygen):

cd doc

doxygen doxyschnaps

firefox html/index.html 



*/












#include "global.h"
#include "geometry.h"
#include "field.h"
#include "pic.h"
#include "skyline.h"
#include "linear_solver.h"
#ifdef _WITH_OPENCL
#include "field_cl.h"
#include "clinfo.h"
#endif

#include "maxwell.h"

#endif
