#include "GaussMethod.h"
#include "node_data.h"
#include "Jacobian.h"

using namespace std;

gauss::gauss(int n) : N(n) {
    if (n == 1) {
        points = { 0.0 };
        weights = { 2.0 };
    }
    else if (n == 2) {
        points = { -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0) };
        weights = { 1.0, 1.0 };
    }
    else if (n == 3) {
        points = { -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0) };
        weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    }
    else {
        throw std::invalid_argument("za duzo punktow");
    }
}

double kwadratura_1D(gauss gauss, double (*func)(double)) {
    double result = 0;
    for (int i = 0; i < gauss.N; i++) {
        result += gauss.weights[i] * func(gauss.points[i]);
    }
    return result;
}

double kwadratura_2D(gauss gauss, double (*func)(double, double)) {
    double result = 0;
    for (int i = 0; i < gauss.N; i++) {
        for (int j = 0; j < gauss.N; j++) {
            result += gauss.weights[i] * gauss.weights[j] * func(gauss.points[i], gauss.points[j]);
        }
    }
    return result;
}

double prostokaty_1D(double a, double b, double N, double(*func)(double)) {
    double h = (b - a) / N; //szerokosc na osi przedzialu
    double result = 0;

    for (int i = 0; i < N; i++) {
        double mid = a + (i + 0.5) * h; //srodek kazdego prostokata
        result += h * func(mid);
    }
    return result;
}

double func1(double x) {
    return 5 * (x * x) + 3 * x + 6;
}

double func2(double x, double y) {
    return 5 * pow(x, 2) * pow(y, 2) + 3 * x * y + 6;
}
