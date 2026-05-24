#include "NumericalSchemes.h"
#include <algorithm>  // std::max
#include <cmath>      // std::abs



namespace NumericalSchemes {
    StateVector getRusanov (Cellule Gauche, Cellule Droite, double gamma) {
        StateVector F_interface;

        //Récupération des vitesses du son (gauche et droite) "Ci"
        double  c_g = Gauche.getSoundSpeed(gamma);
        double  c_d = Droite.getSoundSpeed(gamma);

        //Récupération des vitesses  (gauche et droite) "Ui"
        double u_g = Gauche.getFluidVelocity();
        double u_d = Droite.getFluidVelocity();

        // Récupération des flux des cellules (gauche et droite) "Fi"
        StateVector F_g = Gauche.getFlux(gamma);
        StateVector F_d = Droite.getFlux(gamma);

        // Calcul de Smax
        double S_max = std::max(std::abs(u_g) + c_g , std::abs(u_d) + c_d);


        // Calcul du flux à l'interface

        F_interface.rho = 0.5 * (F_d.rho + F_g.rho) - 0.5 * S_max * (Droite.U.rho - Gauche.U.rho);
        F_interface.rho_u = 0.5 * (F_d.rho_u + F_g.rho_u) - 0.5 * S_max * (Droite.U.rho_u - Gauche.U.rho_u);
        F_interface.rho_E = 0.5 * (F_d.rho_E + F_g.rho_E) - 0.5 * S_max * (Droite.U.rho_E - Gauche.U.rho_E);

        return F_interface;
    }

    // L'outil n°2 (La coquille vide pour l'année prochaine)
    StateVector getRoe(Cellule Gauche, Cellule Droite, double gamma) {
        StateVector F_interface;
        // A FAIRE PLUS TARD : Matrice de Roe
        return F_interface;
    }

    // L'outil n°3 (La coquille vide pour projet de recherche HLLC = standard industriel )
    StateVector getHLLC(Cellule Gauche, Cellule Droite, double gamma) {
        StateVector F_interface;
        // A FAIRE PLUS TARD : Calcul des 3 ondes
        return F_interface;
    }
}
