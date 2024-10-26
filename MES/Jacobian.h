#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <vector>
#include "MES.cpp"

struct Jakobian {
    double J[4][4];
    double J1[4][4];
    double detJ;
};

struct UnivElement {
    std::vector<std::vector<double>> dN_dKsi;
    std::vector<std::vector<double>> dN_dEta;

    UnivElement(int npc); //prototyp kontruktora

    void calcDeriv(); //prototyp metody obliczania pochodnych
};

void computeJacobian(const std::vector<node>& nodes, const std::vector<std::vector<double>>& dN_dKsi, const std::vector<std::vector<double>>& dN_dEta, Jakobian& J);

void displayResults(const UnivElement& uElement, const Jakobian& J);

#endif // !JACOBIAN_H
#pragma once