#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

#include "node_data.h"
#include "GaussMethod.h"
#include "Jacobian.h"

using namespace std;

int main() {
    GlobalData data;
    vector<node> nodes;
    vector<element> elements;
    vector<int> boundary_conditions;

    //ifstream file("dane.txt");
    //if (!file) {
    //    cerr << "brak pliku\n";
    //    return 1;
    //}

    //string line;
    //while (getline(file, line)) {
    //    stringstream ss(line);
    //    string label;
    //    ss >> label;

    //    if (label == "SimulationTime") {
    //        ss >> data.SimulationTime;
    //    }
    //    else if (label == "SimulationStepTime") {
    //        ss >> data.SimulationStepTime;
    //    }
    //    else if (label == "Conductivity") {
    //        ss >> data.Conductivity;
    //    }
    //    else if (label == "Alfa") {
    //        ss >> data.alfa;
    //    }
    //    else if (label == "Tot") {
    //        ss >> data.Tot;
    //    }
    //    else if (label == "InitialTemp") {
    //        ss >> data.InitialTemp;
    //    }
    //    else if (label == "Density") {
    //        ss >> data.Density;
    //    }
    //    else if (label == "SpecificHeat") {
    //        ss >> data.SpecificHeat;
    //    }
    //    else if (label == "Nodes") {
    //        ss.ignore(7); // skip "number"
    //        ss >> data.nN;
    //    }
    //    else if (label == "Elements") {
    //        ss.ignore(7); // skip "number"
    //        ss >> data.nE;
    //    }
    //    else if (label == "*Node") {
    //        while (getline(file, line) && !line.empty() && line[0] != '*') {
    //            int id;
    //            double x, y;
    //            stringstream nodeStream(line);
    //            nodeStream >> id;
    //            nodeStream.ignore(1); // pomin przecinek
    //            nodeStream >> x;
    //            nodeStream.ignore(1);
    //            nodeStream >> y;
    //            nodes.push_back({ x, y });
    //        }
    //    }
    //    else if (label == "*Element") {
    //        while (getline(file, line) && !line.empty() && line[0] != '*') {
    //            stringstream elementStream(line);
    //            int id;
    //            elementStream >> id;

    //            element e;
    //            e.ID = id;

    //            int node_id;
    //            char comma;
    //            while (elementStream >> comma >> node_id) {
    //                e.nodes.push_back(node_id);
    //            }
    //            elements.push_back(e);
    //        }
    //    }
    //    else if (label == "*BC") {
    //        while (getline(file, line) && !line.empty() && line[0] != '*') {
    //            stringstream bcStream(line);
    //            int nodeID;
    //            while (bcStream >> nodeID) {
    //                boundary_conditions.push_back(nodeID);
    //                if (bcStream.peek() == ',') bcStream.ignore();
    //            }
    //        }
    //    }
    //}
    //file.close();

    //cout << "Global Data:" << endl;
    //cout << "SimulationTime: " << data.SimulationTime << endl;
    //cout << "SimulationStepTime: " << data.SimulationStepTime << endl;
    //cout << "Conductivity: " << data.Conductivity << endl;
    //cout << "Alfa: " << data.alfa << endl;
    //cout << "Tot: " << data.Tot << endl;
    //cout << "InitialTemp: " << data.InitialTemp << endl;
    //cout << "Density: " << data.Density << endl;
    //cout << "SpecificHeat: " << data.SpecificHeat << endl;
    //cout << "Number of Nodes: " << data.nN << endl;
    //cout << "Number of Elements: " << data.nE << endl;

    //cout << "\nNodes:" << endl;
    //for (const auto& n : nodes) {
    //    cout << "x: " << n.x << ", y: " << n.y << endl;
    //}

    //cout << "\nElements:" << endl;
    //for (const auto& e : elements) {
    //    cout << "Element ID: " << e.ID << " Nodes: ";
    //    for (int id : e.nodes) cout << id << " ";
    //    cout << endl;
    //}

    //cout << "\nBoundary Conditions:" << endl;
    //for (int bc : boundary_conditions) cout << bc << " ";
    //cout << endl;

    //cout << fixed << setprecision(10);

    //cout << "kwadratura dla 1D\n";
    //for (int i = 1; i <= 3; i++) {
    //    gauss gauss1D(i);
    //    double result = kwadratura_1D(gauss1D, func1);
    //    cout << "Dla N = " << i << ": wynik = " << result << endl;
    //}

    //cout << "kwadratura dla 2D\n";
    //for (int i = 1; i <= 3; i++) {
    //    gauss gauss2D(i);
    //    double result = kwadratura_2D(gauss2D, func2);
    //    cout << "Dla N = " << i << ": wynik = " << result << endl;
    //}

    //cout << "metoda prostokatow \n";
    //cout << prostokaty_1D(-1, 1, 10, func1) << endl;


    vector<node> elementNodes = { {0, 0}, {0.025, 0}, {0.025, 0.025}, {0, 0.025} };

    
    UnivElement uElement(4);
    uElement.calcDeriv();


    displayResults(uElement, elementNodes);

    return 0;
}