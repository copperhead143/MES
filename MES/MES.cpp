#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

#include "GaussMethod.h"
#include "Jacobian.h"


using namespace std;

struct node {
    double x;
    double y;
};

struct element {
    int ID;
    vector<int> nodes;
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



int main() {
    GlobalData data;
    vector<node> nodes;
    vector<element> elements;
    vector<int> boundary_conditions;

    ifstream file("dane.txt");
    if (!file) {
        cerr << "brak pliku\n";
        return 1;
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string label;
        ss >> label;

        if (label == "SimulationTime") {
            ss >> data.SimulationTime;
            //cout << "Loaded SimulationTime: " << data.SimulationTime << endl;
        }
        else if (label == "SimulationStepTime") {
            ss >> data.SimulationStepTime;
            //cout << "Loaded SimulationStepTime: " << data.SimulationStepTime << endl;
        }
        else if (label == "Conductivity") {
            ss >> data.Conductivity;
            //cout << "Loaded Conductivity: " << data.Conductivity << endl;
        }
        else if (label == "Alfa") {
            ss >> data.alfa;
            //cout << "Loaded Alfa: " << data.alfa << endl;
        }
        else if (label == "Tot") {
            ss >> data.Tot;
            //cout << "Loaded Tot: " << data.Tot << endl;
        }
        else if (label == "InitialTemp") {
            ss >> data.InitialTemp;
            //cout << "Loaded InitialTemp: " << data.InitialTemp << endl;
        }
        else if (label == "Density") {
            ss >> data.Density;
            //cout << "Loaded Density: " << data.Density << endl;
        }
        else if (label == "SpecificHeat") {
            ss >> data.SpecificHeat;
            //cout << "Loaded SpecificHeat: " << data.SpecificHeat << endl;
        }
        else if (label == "Nodes") {
            ss.ignore(7); // skip "number"
            ss >> data.nN;
            //cout << "Loaded Number of Nodes: " << data.nN << endl;
        }
        else if (label == "Elements") {
            ss.ignore(7); // skip "number"
            ss >> data.nE;
            //cout << "Loaded Number of Elements: " << data.nE << endl;
        }
        else if (label == "*Node") {
            //cout << "Entering Node section" << endl;
            //odczyt node
            while (getline(file, line) && !line.empty() && line[0] != '*') { //dopoki jest linia, nie jest pusta i pierwszy znak to nie *
                int id;
                double x, y;
                stringstream nodeStream(line);
                nodeStream >> id;
                nodeStream.ignore(1); // pomin przcinek
                nodeStream >> x;
                nodeStream.ignore(1);
                nodeStream >> y;
                nodes.push_back({ x, y });
                //cout << "Loaded Node - ID: " << id << ", x: " << x << ", y: " << y << endl;
            }
        }

        else if (label == "*Element") {

            while (getline(file, line) && !line.empty() && line[0] != '*') {
                stringstream elementStream(line);
                int id;
                elementStream >> id;  //odczyt IDA

                element e;
                e.ID = id;

                int node_id;
                char comma;
                while (elementStream >> comma >> node_id) {
                    e.nodes.push_back(node_id);
                }

                elements.push_back(e);
                //cout << "Loaded Element - ID: " << e.ID << " with nodes: ";
                for (int nodeID : e.nodes) cout << nodeID << " ";
                cout << endl;
            }
        }


        else if (label == "*BC") {
            cout << "Entering Boundary Conditions section" << endl;
            while (getline(file, line) && !line.empty() && line[0] != '*') {
                stringstream bcStream(line);
                int nodeID;
                while (bcStream >> nodeID) {
                    boundary_conditions.push_back(nodeID);
                    if (bcStream.peek() == ',') bcStream.ignore();
                    //cout << "Loaded Boundary Condition Node ID: " << nodeID << endl;
                }
            }
        }
    }
    file.close();

    cout << "Global Data:" << endl;
    cout << "SimulationTime: " << data.SimulationTime << endl;
    cout << "SimulationStepTime: " << data.SimulationStepTime << endl;
    cout << "Conductivity: " << data.Conductivity << endl;
    cout << "Alfa: " << data.alfa << endl;
    cout << "Tot: " << data.Tot << endl;
    cout << "InitialTemp: " << data.InitialTemp << endl;
    cout << "Density: " << data.Density << endl;
    cout << "SpecificHeat: " << data.SpecificHeat << endl;
    cout << "Number of Nodes: " << data.nN << endl;
    cout << "Number of Elements: " << data.nE << endl;

    cout << "\nNodes:" << endl;
    for (const auto& n : nodes) {
        cout << "x: " << n.x << ", y: " << n.y << endl;
    }

    cout << "\nElements:" << endl;
    for (const auto& e : elements) {
        cout << "Element ID: " << e.ID << " Nodes: ";
        for (int id : e.nodes) cout << id << " ";
        cout << endl;
    }

    cout << "\nBoundary Conditions:" << endl;
    for (int bc : boundary_conditions) cout << bc << " ";
    cout << endl;

    //=============================================================================//
    //obliczanie kwadratury gaussa dla dwoch funkcji

    cout << fixed << setprecision(10);

    cout << "kwadratura dla 1D\n";
    for (int i = 1; i <= 3; i++) {
        gauss gauss1D(i);
        double result = kwadratura_1D(gauss1D, func1);
        cout << "Dla N = " << i << ": wynik = " << result << endl;
    }


    cout << "kwadratura dla 2D\n";
    for (int i = 1; i <= 3; i++) {
        gauss gauss2D(i);
        double result = kwadratura_2D(gauss2D, func2);
        cout << "Dla N = " << i << ": wynik = " << result << endl;
    }


    cout << "metoda prostokatow \n";
    cout << prostokaty_1D(-1, 1, 10, func1);

    vector<node> nodes1;
    nodes1[0].x = 0;
    nodes1[0].y = 0;
    nodes1[1].x = 0.025;
    nodes1[1].y = 0.0;
    nodes1[2].x = 0.025;
    nodes1[2].y = 0.025;
    nodes1[3].x = 0.0;
    nodes1[3].y = 0.025;



    return 0;
}
