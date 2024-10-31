#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <vector>
#include "node_data.h"

struct Jakobian {
    double J[2][2];
    double J1[2][2];
    double detJ;
};

struct UnivElement {
    std::vector<std::vector<double>> dN_dKsi;
    std::vector<std::vector<double>> dN_dEta;

    UnivElement(int npc);
    void calcDeriv();
};


void computeJacobian(const std::vector<node>& nodes,
    const std::vector<std::vector<double>>& dN_dKsi,
    const std::vector<std::vector<double>>& dN_dEta,
    Jakobian& J, int gaussPointIndex);

void computeDN_dx_dy(const UnivElement& uElement, const Jakobian& J, std::vector<std::vector<double>>& dN_dx, std::vector<std::vector<double>>& dN_dy);

void displayResults(const UnivElement& uElement, const std::vector<node>& elementNodes);

#endif