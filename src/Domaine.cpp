
#include "Domaine.h"

#include <iostream>
#include <fstream>
#include <limits>    // std::numeric_limits
#include <algorithm> // std::min
#include <cmath>     // std::abs

double Domaine::getDeltat (const Cellule& cell, double gamma, double CFL) const{
    // initialisation de c (la vietsse du son)
    const double c = cell.getSoundSpeed(gamma) ;

    // récupération de la vitesse
    const double u = cell.getFluidVelocity() ;

    // Calcul du delta t
    double delta_t = CFL * cell.dx / (std::abs(u) + c); ;

    return delta_t ;
}

double Domaine::getGlobalDeltat(std::vector<Cellule> tube, double CFL, double gamma) const {
    double dt_min = 1000000.0;
    int N = tube.size();

    for (int i = 0; i < N; i++) {
        double dt_local = getDeltat(tube[i], gamma, CFL);
        if (dt_local < dt_min) {
            dt_min = dt_local;
        }
    }
    return dt_min;
}

void Domaine::exporterInstantane(int numero_photo) {
    // Crée un fichier nommé resultat_0.csv, resultat_1.csv, etc.
    std::ofstream fichier("resultat_" + std::to_string(numero_photo) + ".csv");

    // En-tête du fichier CSV
    fichier << "x,rho,vitesse,pression\n";

    int N = tube.size() - 2;
    for (int i = 1; i <= N; i++) {
        double u = tube[i].getFluidVelocity();
        // Pression = (rho_E - 0.5 * rho * u^2) * (gamma - 1)
        double p = tube[i].getPressure(gamma);

        // On écrit les colonnes séparées par des virgules
        fichier << tube[i].x << ","
                << tube[i].U.rho << ","
                << u << ","
                << p << "\n";
    }
    fichier.close();
}

// Conditions aux limites transmissives (zéro gradient)
// Cellule fantôme = copie de la cellule réelle voisine
void Domaine::appliquerConditionsLimites() {

    int N = tube.size() - 2;
    tube[0] = tube[1];
    tube[N+1] = tube[N];
}

//Condition Initiale
Domaine::Domaine(int N, double L, double gamma, double t_max, double r) {
    this -> gamma = gamma ;
    this -> t_max = t_max ;

    tube.resize(N+2); // on prend en compte les cellules "fantomes" pour le schema d'ou N+2
    double dx = L/N ;
    for (int i = 1; i < N+1; i++) {
        tube[i].dx = dx ;
        tube[i].x = (i - 0.5) * dx ;

        if (tube[i].x > 1.0) {
            // DROITE : Air froid (Basse pression)
            double T_droite = 300.0;
            tube[i].U.rho = 0.125;
            tube[i].U.rho_u = 0.0;
            double P_droite = tube[i].U.rho * r * T_droite;
            tube[i].U.rho_E = P_droite / (gamma - 1.0);
        }
        else {
            // GAUCHE : Air chaud (Haute pression)
            double T_gauche = 1000.0;
            tube[i].U.rho = 1.0;
            tube[i].U.rho_u = 0.0;
            double P_gauche = tube[i].U.rho * r * T_gauche;
            tube[i].U.rho_E = P_gauche / (gamma - 1.0);
        }

    }
    // Application des conditions aux limites
    appliquerConditionsLimites();

}


void Domaine::run() {
    int N = tube.size() - 2; // Le nombre de vraies cellules
    std::vector<StateVector> flux_interfaces(N + 1);
    double t = 0.0;

    //pour le fichier
    int compteur_iterations = 0;
    int compteur_photos = 0;

    // Photo initiale à t = 0
    exporterInstantane(compteur_photos);
    compteur_photos++;

    while (t < t_max) {
        // Recalcul du raccord de sécurité des fantômes à chaque itération
        appliquerConditionsLimites();
        // Le chronomètre est calculé une seule fois par tour
        double dt = getGlobalDeltat(tube, 0.5, gamma);

        // Étape 1 : La photographie des flux aux interfaces
        for (int i = 0; i <= N; i++) {
            flux_interfaces[i] = NumericalSchemes::getRusanov(tube[i], tube[i+1], gamma);
        }

        // Étape 2 : La mise à jour des vraies usines (de 1 à N)
        for (int i = 1; i <= N; i++) {
            tube[i].U.rho   = tube[i].U.rho   - (dt/tube[i].dx) * (flux_interfaces[i].rho   - flux_interfaces[i-1].rho);
            tube[i].U.rho_u = tube[i].U.rho_u - (dt/tube[i].dx) * (flux_interfaces[i].rho_u - flux_interfaces[i-1].rho_u);
            tube[i].U.rho_E = tube[i].U.rho_E - (dt/tube[i].dx) * (flux_interfaces[i].rho_E - flux_interfaces[i-1].rho_E);
        }

        // Le temps avance !
        t = t + dt;

        compteur_iterations++;

        // Toutes les 5 itérations, on prend une photo !
        if (compteur_iterations % 5 == 0) {
            exporterInstantane(compteur_photos);
            compteur_photos++;
        }
    }
    std::cout << "Simulation terminee avec succes a t = " << t << " s !" << std::endl;
    // Une dernière photo à la fin
    exporterInstantane(compteur_photos);
    std::cout << "Simulation terminee ! " << compteur_photos << " fichiers CSV generes !" << std::endl;
}