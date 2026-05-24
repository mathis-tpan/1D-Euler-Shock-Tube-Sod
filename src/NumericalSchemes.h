#pragma once
#include "Cellule.h"


namespace NumericalSchemes {

    // L'outil n°1 (Rusanov):
    StateVector getRusanov (Cellule Gauche, Cellule Droite, double gamma) ;

    // L'outil n°2 (La coquille vide pour l'année prochaine)
    StateVector getRoe(Cellule Gauche, Cellule Droite, double gamma);

    // L'outil n°3 (La coquille vide pour projet de recherche HLLC = standard industriel )
    StateVector getHLLC(Cellule Gauche, Cellule Droite, double gamma) ;

};