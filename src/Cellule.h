#pragma once
#include "StateVector.h"

struct Cellule {
    double dx ;
    double x ;
    StateVector U ;

    double getPressure(double gamma) const ;
    double getTemperature( double gamma, double r) const ;
    double getSoundSpeed(double gamma) const ;
    double getFluidVelocity() const ;
    StateVector getFlux( double gamma);
};