#pragma once
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#define h (pow(10, -5))
#define g (1.0-1.0/sqrt(2))
#define b (1.0-1.0/sqrt(2))
#define p1 (1.0-1.0/sqrt(2))
#define p2 (1.0/sqrt(2))
#define L_4 pow(0.05, 4)
#define G (6.582 / (1.22 * 1.22) * pow(10, -51))
#define a0  (2 * sqrt(24 * 3.14159 * G) * pow(10, 23))

int main()
{
    double a = 5.0;
    double f = pow(10, 23);
    //double H_inside = (2 * a * a * a - a) * pow(L_4 * z * z / 2.0, a * a) + L_4 * (1 - cos(y));
    //if (H_inside < 0) H_inside = 0;  // Защита от NaN
    //double H = sqrt(H_inside);
    //double dzdt = -sqrt(24 * M_PI * G) * f * z / (2 * a * a - 1) * (H / 0.0025) - 1 / (2 * pow(a, 5) - a * a * a) * sin(y) * pow(L_4 * z * z / 2.0, 1 - a * a);
    cout << "A= " << sqrt(24 * M_PI * G) * f / (2 * a * a - 1)/0.0025 << endl;
    cout << "B= " << (2 * a * a * a - a) * pow (L_4/2.0, a*a) << endl;
    double B = (2 * a * a * a - a) * pow(L_4 / 2.0, a * a);
    cout << "C= " << L_4 << endl;
    cout << "D= " << 1 / (2 * pow(a, 5) - a * a * a) * pow(L_4 / 2.0, 1 - a * a) << endl << endl;
    double D = 1 / (2 * pow(a, 5) - a * a * a) * pow(L_4 / 2.0, 1 - a * a);
    cout << "z= " << pow(10, 128.0 / 48.0) << "*x" << endl << endl;
    double k = pow(10, 128.0 / 48.0);
    cout << "A= " << sqrt(24 * M_PI * G) * f / (2 * a * a - 1) / 0.0025 << endl;
    cout << "B'= " << pow(k, 50) * B << endl;
    cout << "C= " << L_4 << endl;
    cout << "D'= " << D/pow(k, 50)  << endl;
}
