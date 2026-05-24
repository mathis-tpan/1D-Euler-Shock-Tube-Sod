#include "src/Domaine.h"

int main() {
    // 1. Paramètres de la simulation
    int N = 500;
    double L = 2.0;
    double gamma = 1.4;
    double r = 287.0;
    double t_max = 0.002;

    // définir le domaine 1D dans notre cas !
    Domaine monTube(N, L, gamma, t_max, r);

    // on lance la simulation
    monTube.run();

    return 0;
}