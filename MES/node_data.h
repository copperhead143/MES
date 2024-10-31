#pragma once
#include <vector>

struct node {
    double x;
    double y;
};

struct element {
    int ID;
    std::vector<int> nodes;
};

struct GlobalData {
    double SimulationTime;
    double SimulationStepTime;
    double Conductivity;
    double alfa;
    double Tot;
    double InitialTemp;
    double Density;
    double SpecificHeat;
    int nN; // number of nodes
    int nE; // number of elements
    double npc; //ilosc punktow calkowania
};

