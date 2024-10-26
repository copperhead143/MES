#ifndef GAUSSMETHOD_H
#define GAUSSMETHOD_H

#include <vector>
#include <stdexcept>
#include <cmath>

struct gauss {
    std::vector<double> points;  
    std::vector<double> weights; 
    int N;                       

    gauss(int n); //prototyp konstruktora
};

//prototypy funkcji
double kwadratura_1D(gauss gauss, double (*func)(double));
double kwadratura_2D(gauss gauss, double (*func)(double, double));
double prostokaty_1D(double a, double b, double N, double(*func)(double));
double func1(double x);
double func2(double x, double y);

#endif
#pragma once
