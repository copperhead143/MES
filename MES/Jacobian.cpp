#include "Jacobian.h"
#include "MES.cpp"

#include <vector>

UnivElement::UnivElement(int npc) {
    dN_dKsi.resize(npc, std::vector<double>(4));
    dN_dEta.resize(npc, std::vector<double>(4));
}

void UnivElement::calcDeriv() {
    for (int i = 0; i < dN_dKsi.size(); ++i) {
        double ksi = dN_dKsi[i][0];
        double eta = dN_dEta[i][0];

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

void computeJacobian(const std::vector<node>& nodes, const std::vector<std::vector<double>>& dN_dKsi, const std::vector<std::vector<double>>& dN_dEta, Jakobian& J) {
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            J.J[i][j] = 0;
            for (int k = 0; k < 4; ++k) {
                J.J[0][0] += dN_dKsi[k][i] * nodes[k].x;
                J.J[0][1] += dN_dKsi[k][i] * nodes[k].y;
                J.J[1][0] += dN_dEta[k][j] * nodes[k].x;
                J.J[1][1] += dN_dEta[k][j] * nodes[k].y;
            }
        }
    }
    //liczenie wyznacznika
    J.detJ = J.J[0][0] * J.J[1][1] - J.J[0][1] * J.J[1][0];

    //odwrócenie jakobianu
    J.J1[0][0] = J.J[1][1] / J.detJ;
    J.J1[0][1] = -J.J[0][1] / J.detJ;
    J.J1[1][0] = -J.J[1][0] / J.detJ;
    J.J1[1][1] = J.J[0][0] / J.detJ;

}


void displayResults(const UnivElement& uElement, const Jakobian& J) {
    cout << fixed << setprecision(5);
    cout << "Tabela pochodnych funkcji kszta³tu dla punktów ca³kowania Gaussa:\n";
    cout << "dN/dKsi\t\t\tdN/dEta\n";
    for (int i = 0; i < 4; ++i) {
        cout << "pc" << (i + 1) << " ";
        for (int j = 0; j < 4; ++j) {
            cout << uElement.dN_dKsi[i][j] << "\t" << uElement.dN_dEta[i][j] << "\t";
        }
        cout << endl;
    }

    cout << "\nMacierz Jakobianu:\n";
    cout << "J = [" << J.J[0][0] << ", " << J.J[0][1] << "; " << J.J[1][0] << ", " << J.J[1][1] << "]\n";
    cout << "detJ = " << J.detJ << endl;

    cout << "\nOdwrócony Jakobian:\n";
    cout << "J1 = [" << J.J1[0][0] << ", " << J.J1[0][1] << "; " << J.J1[1][0] << ", " << J.J1[1][1] << "]\n";
}