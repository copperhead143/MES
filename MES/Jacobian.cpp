#include <vector>
#include <iomanip>
#include <iostream>

#include "Jacobian.h"
#include "node_data.h"
#include "GaussMethod.h"


using namespace std;

UnivElement::UnivElement(int npc) {
    dN_dKsi.resize(npc, std::vector<double>(4));
    dN_dEta.resize(npc, std::vector<double>(4));
}

void UnivElement::calcDeriv() {
    std::vector<std::pair<double, double>> gauss_points = {
        {-1.0 / sqrt(3.0), -1.0 / sqrt(3.0)},
        {1.0 / sqrt(3.0), -1.0 / sqrt(3.0)},
        {1.0 / sqrt(3.0), 1.0 / sqrt(3.0)},
        {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}
    };

    dN_dKsi.resize(4, std::vector<double>(4));
    dN_dEta.resize(4, std::vector<double>(4));

    for (int i = 0; i < 4; ++i) {
        double ksi = gauss_points[i].first;
        double eta = gauss_points[i].second;

        dN_dKsi[i][0] = -0.25 * (1 - eta);
        dN_dKsi[i][1] = 0.25 * (1 - eta);
        dN_dKsi[i][2] = 0.25 * (1 + eta);
        dN_dKsi[i][3] = -0.25 * (1 + eta);

        dN_dEta[i][0] = -0.25 * (1 - ksi);
        dN_dEta[i][1] = -0.25 * (1 + ksi);
        dN_dEta[i][2] = 0.25 * (1 + ksi);
        dN_dEta[i][3] = 0.25 * (1 - ksi);
    }
}

void computeJacobian(const std::vector<node>& nodes, const std::vector<std::vector<double>>& dN_dKsi, const std::vector<std::vector<double>>& dN_dEta, Jakobian& J, int gaussPointIndex) {
    J.J[0][0] = J.J[0][1] = J.J[1][0] = J.J[1][1] = 0;

    for (int i = 0; i < 4; ++i) {
        J.J[0][0] += dN_dKsi[gaussPointIndex][i] * nodes[i].x;
        J.J[0][1] += dN_dKsi[gaussPointIndex][i] * nodes[i].y;
        J.J[1][0] += dN_dEta[gaussPointIndex][i] * nodes[i].x;
        J.J[1][1] += dN_dEta[gaussPointIndex][i] * nodes[i].y;
    }

    J.detJ = J.J[0][0] * J.J[1][1] - J.J[0][1] * J.J[1][0];

    J.J1[0][0] = J.J[1][1] / J.detJ;
    J.J1[0][1] = -J.J[0][1] / J.detJ;
    J.J1[1][0] = -J.J[1][0] / J.detJ;
    J.J1[1][1] = J.J[0][0] / J.detJ;
}


void computeDN_dx_dy(const UnivElement& uElement, const Jakobian& J, std::vector<std::vector<double>>& dN_dx, std::vector<std::vector<double>>& dN_dy) {
    dN_dx.resize(4, std::vector<double>(4));
    dN_dy.resize(4, std::vector<double>(4));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) { 
            dN_dx[i][j] = J.J1[0][0] * uElement.dN_dKsi[i][j] + J.J1[0][1] * uElement.dN_dEta[i][j];
            dN_dy[i][j] = J.J1[1][0] * uElement.dN_dKsi[i][j] + J.J1[1][1] * uElement.dN_dEta[i][j];
        }
    }
}


void displayResults(const UnivElement& uElement, const std::vector<node>& elementNodes) {
    cout << fixed << setprecision(5);

    for (int i = 0; i < 4; ++i) {
        Jakobian J;
        computeJacobian(elementNodes, uElement.dN_dKsi, uElement.dN_dEta, J, i);

        std::vector<std::vector<double>> dN_dx, dN_dy;
        computeDN_dx_dy(uElement, J, dN_dx, dN_dy);

        cout << "=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=\n";
        cout << "Punkt Gaussa " << (i + 1) << ":\n";
        cout << "dN/dKsi\t\tdN/dEta\n";
        for (int j = 0; j < 4; ++j) {
            cout << uElement.dN_dKsi[i][j] << "\t" << uElement.dN_dEta[i][j] << "\n";
        }

        cout << "\nMacierz Jacobiego:\n";
        cout << "J = [" << J.J[0][0] << ", " << J.J[0][1] << "; " << J.J[1][0] << ", " << J.J[1][1] << "]\n";
        cout << "detJ = " << J.detJ << endl;

        cout << "\nOdwrocony Jacobian:\n";
        cout << "J1 = [" << J.J1[0][0] << ", " << J.J1[0][1] << "; " << J.J1[1][0] << ", " << J.J1[1][1] << "]\n\n";

        cout << "Pochodne wzglêdem x i y:\n";
        for (int j = 0; j < 4; ++j) {
            cout << "dN" << (j + 1) << "/dx: " << dN_dx[i][j] << "\t";
            cout << "dN" << (j + 1) << "/dy: " << dN_dy[i][j] << "\n";
        }
        cout << endl;
    }
}