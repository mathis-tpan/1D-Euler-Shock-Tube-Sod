#include "Cellule.h"
#include <cmath>

#include "StateVector.h"

double Cellule::getPressure(double gamma) const {
    double p ;
    p = (gamma - 1.0) * (U.rho_E  - (U.rho_u * U.rho_u)/(2.0*U.rho));
    return p;
}

double Cellule::getTemperature( double gamma, double r) const  {
    double p = getPressure(gamma);

    double T = p / (U.rho * r);
    return T;
}

double Cellule::getSoundSpeed(double gamma) const {
    double p = getPressure(gamma);
    return std::sqrt((gamma * p) / U.rho);
}

double Cellule::getFluidVelocity() const {
    return U.rho_u / U.rho;
}

StateVector Cellule::getFlux( double gamma) {
    StateVector F;

    double p = getPressure(gamma);
    double u = getFluidVelocity();

    //calcul des flux
    F.rho = U.rho_u ; // flux de masse
    F.rho_u = U.rho_u * u + p  ; // flux de quantité de mouvement
    F.rho_E = u * (U.rho_E + p) ; // flux d'energie

    return F; // on retourne le vecteur flux F
}