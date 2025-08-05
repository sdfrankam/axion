#pragma once
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#define h (pow(10, -11))
#define g (1.0-1.0/sqrt(2))
#define b (1.0-1.0/sqrt(2))
#define p1 (1.0-1.0/sqrt(2))
#define p2 (1.0/sqrt(2))
#define L_4 pow(0.05, 4)
#define G (6.582 / (1.22 * 1.22) * pow(10, -51))
#define a0  (2 * sqrt(24 * 3.14159 * G) * pow(10, 23))

#define k 464.159
#define A 0.471372
#define B 0.0124094
#define C 6.25e-06
#define D 1.0073e-05

//void system_ODE_q(double z, double y) {
//    double dzdt, dydt;
//    dzdt = f_q(z, y);
//    dydt = z;
//}

double f_test(double z, double y) {
    return -y;
}

void Jacoby_Matrics_test(vector<double>& J, double z, double y) {
    //   | J[0] J[1] |
    // J=|           |
    //   | J[2] J[3] |
    J[0] = 0;
    J[1] = -1;
    J[2] = 1.0;
    J[3] = 0;
}

double f_q(double z, double y) {
    //G=hc/(m_p)^2=] c=1 [ =h/(m_p)^2 = 
    double a = sqrt(8 * 3.14159 * G * 3) * pow(10, 23);
    return -sin(y) - a * z * sqrt(z * z / 2.0 - cos(y) + 1);
}

void Jacoby_Matrics_q(vector<double>& J, double z, double y) {
    //   | J[0] J[1] |
    // J=|           |
    // double a = sqrt(8 * 3.14159 * G * 3) * pow(10, 23);
    //   | J[2] J[3] |
    double a = sqrt(8 * 3.14159 * G * 3) * pow(10, 23);
    J[0] = -a * sqrt(z * z / 2.0 - cos(y) + 1) - a * z * z / (2.0 * sqrt(z * z / 2.0 - cos(y) + 1));
    J[1] = -cos(y) - a * z * sin(y) / (2 * sqrt(z * z / 2.0 - cos(y) + 1));
    J[2] = 1.0;
    J[3] = 0;
}


//void system_ODE_f(double z, double y) {
//    double dzdt, dydt;
//    dzdt = f_f(z, y);
//    dydt = z;
//}

double f_f(double z, double y) {
    //G=hc/(m_p)^2=] c=1 [ =h/(m_p)^2 = 
    double a = sqrt(8 * 3.14159 * G * 3) * pow(10, 23);
    return sin(y) - a * z * sqrt(-z * z / 2.0 - cos(y) + 1);
}

static void Jacoby_Matrics_f(vector<double>& J, double z, double y) {
    //   | J[0] J[1] |
    // J=|           |
    //   | J[2] J[3] |
    double a = sqrt(8 * 3.14159 * G * 3) * pow(10, 23);
    J[0] = -a * sqrt(-z * z / 2.0 - cos(y) + 1) + a * z * z / (2.0 * sqrt(-z * z / 2.0 - cos(y) + 1));
    J[1] = cos(y) - a * z * sin(y) / (2 * sqrt(-z * z / 2.0 - cos(y) + 1));
    J[2] = 1.0;
    J[3] = 0;
}

double f_t(double z, double y) {
    //G=hc/(m_p)^2=] c=1 [ =h/(m_p)^2 = 
    double a = sqrt(8 * 3.14159 * G * 3) * pow(10, 23);
    return (L_4 * z * z - 1) * (a * z * sqrt((1 - cos(y)) / 2 * (2 - L_4 * z * z) / (1 - L_4 * z * z)) + 1 / (tan(y / 2) * L_4));
}
static void Jacoby_Matrics_t(vector<double>& J, double z, double y) {
    //   | J[0] J[1] |
    // J=|           |

    double a = sqrt(8 * 3.14159 * G * 3) * pow(10, 23);

    double E1 = (2 - L_4 * z * z) / ((1 - L_4 * z * z) * sqrt(1 - L_4 * z * z)) - 1 / sqrt((2 - L_4 * z * z) * (1 - L_4 * z * z));
    double E2 = sqrt((1 - cos(y)) / 2 * (2 - L_4 * z * z) / (1 - L_4 * z * z)) + 
        L_4 * z * z * sqrt((1 - cos(y)) / 2) * E1;
    double E3 = 2 * L_4 * z * 
        (a * z * sqrt((1 - cos(y)) / 2 * (2 - L_4 * z * z) / (1 - L_4 * z * z))+1/(tan(y/2)*L_4)) 
        + (L_4 * z * z - 1) * a * E2;

    J[0] = E3;

    J[1] = (L_4*z*z-1)*(a*z*sqrt((2-L_4*z*z)/(1-L_4*z*z))-1/(2*sin(y/2)*sin(y/2)*L_4));
        
    J[2] = 1.0;
    J[3] = 0;
}

//static void Jacoby_Matrics_t(vector<double>& J, double z, double y) {
//    //   | J[0] J[1] |
//    // J=|           |
//
//    double a = sqrt(8 * 3.14159 * G * 3) * pow(10, 23);
//
//    double E1 = (2 - L_4 * z * z) / ((1 - L_4 * z * z) * sqrt(1 - L_4 * z * z)) - 1 / sqrt((2 - L_4 * z * z) * (1 - L_4 * z * z));
//    double E2 = sqrt((1 - cos(y)) / 2 * (2 - L_4 * z * z) / (1 - L_4 * z * z)) +
//        L_4 * z * z * sqrt((1 - cos(y)) / 2) * E1;
//    double E3 = 2 * L_4 * z *
//        (a * z * sqrt((1 - cos(y)) / 2 * (2 - L_4 * z * z) / (1 - L_4 * z * z)) + 1 / (tan(y / 2) * L_4))
//        + (L_4 * z * z - 1) * a * E2;
//
//    J[0] = E3;
//
//    J[1] = (L_4 * z * z - 1) * (a * z * sqrt((2 - L_4 * z * z) / (1 - L_4 * z * z)) - 1 / (2 * sin(y / 2) * sin(y / 2) * L_4));
//
//    J[2] = 1.0;
//    J[3] = 0;
//}

//double f_generiliezed_1(double z, double y) {
//    double a=10.0, f=pow(10, 23);
//    double X = L_4 * z * z / 2.0;
//    double dzdt = -(sqrt(8 * M_PI * G) * pow(a, 3) * f * pow(X, pow(a, 2)) * sqrt((2 * pow(a, 3)-a) * pow(X, pow(a, 2)) + L_4 * (1.0 - cos(y))) + L_4 * 0.05 *0.05 * z * sin(y)) / (L_4 * (2 * pow(a, 5) - pow(a, 3)) * pow(X, pow(a, 2)));
//    return dzdt;
//}
//
//void Jacoby_Matrix_f_g_1(vector<double>& J, double z, double y) {
//    double a = 1.0, f = pow(10, 23);
//    double X= L_4 * z * z / 2.0;
//    double A = sqrt(24 * M_PI * G * 2 * pow(a, 3) * f * pow(L_4 / 2, pow(a, 2)));
//    double B = sqrt((2 * pow(a, 3) - a) * pow(L_4 * z * z / 2, pow(a, 2)) + L_4 * pow(1 - cos(y), 2));
//    double C = (2 * pow(a, 3) - a) * pow(L_4 / 2, pow(a, 2));
//    double D = pow(L_4, a * a + 1) * (2 * pow(a, 5) - pow(a, 3)) / pow(2, a * a);
//    double G1 = (sqrt(8 * M_PI * G) * pow(a, 3) * f * pow(X, pow(a, 2)) * sqrt((2 * pow(a, 3) - a) * pow(X, pow(a, 2)) + L_4 * (1.0 - cos(y))) + L_4 * 0.05*0.05 * z * sin(y));
//    double G2 = (L_4 * (2 * pow(a, 5) - pow(a, 3)) * pow(X, pow(a, 2)));
//
//    J[0] = -(A * 2 * pow(a, 2) * pow(z, 2 * pow(a, 2) - 1) * B + A * pow(z, 2 * pow(a, 2) - 1) * (C * pow(z, 2 * pow(a, 2) - 1) * 2 * pow(a, 2)) / (2 * B) + pow(L_4, 2) * sin(y)) / (L_4 * (2 * pow(a, 5) - pow(a, 3)) * pow(L_4 * z * z / 2, pow(a, 2)))+2*pow(a, 2)* D *pow(z, a*a-1)*G1/(G2*G2);
//
//    double A1 = sqrt(24 * M_PI * G / 3) * 2 * pow(a, 3) + pow(L_4 * z * z / 2, pow(a, 2));
//    double B1 = (2 * pow(a, 3) - a) * pow(L_4 * z * z / 2, pow(a, 2));
//    double D1 = L_4 * (2 * pow(a, 5) - pow(a, 3)) * pow(L_4 * z * z / 2, pow(a, 2));
//
//    J[1] = -1.0 / D1 * (L_4 * sin(y) * A1 / sqrt(B1 + L_4 * (1 - cos(y))) + pow(0.05, 6)* cos(y)); //изменила pow(L_4, 8) при косинусе на 0.05^6
//    J[2] = 1.0;
//    J[3] = 0;
//
//}
//double g_1(double z, double y) {
//    double a = 3.0;
//    double f = pow(10, 23);    
//    double H_inside = (2 * a * a * a - a) * pow(L_4 * z * z / 2.0, a * a) + L_4 * (1 - cos(y));
//    if (H_inside < 0) H_inside = 0;  // Защита от NaN
//    double H = sqrt(H_inside);
//    double dzdt = -sqrt(24 * M_PI * G) * f * z / (2 * a * a - 1) * (H/0.0025) - 1 / (2 * pow(a, 5) - a * a * a) * sin(y) * pow(L_4 * z * z / 2.0, 1 - a * a);
//    return dzdt;
//}
//void Jacoby_Matrics_g_1(vector<double>& J, double z, double y) {
//    double a = 3.0, f = pow(10, 23);
//
//    double H = sqrt((2 * a * a * a - a) * pow(L_4 * z * z / 2.0, a * a) + L_4 * (1 - cos(y)));
//
//    J[0] = -f * sqrt(24 * M_PI * G) / ((2 * a * a - 1) * 0.05 * 0.05) * H - a * a * a * pow(L_4 * z * z / 2.0, a * a) * f * sqrt(24 * M_PI * G) / (0.05 * 0.05 * H) - (2 - 2 * a * a) / (2 * pow(a, 5) - a * a * a) * sin(y) * pow(L_4 / 2.0, 1 - a * a) */* это слагаемое и срет гадостью (z*z)^(1/2-1) */pow(z, 1 - a * a * 2.0);
//    
//    J[1] = -z * sqrt(24 * M_PI * G) * f / (2 * a * a - 1) * 0.05 * 0.05 * sin(y) / (2 * H) - 1 / (2 * pow(a, 5) - a * a * a) * cos(y) * pow(L_4 * z * z / 2.0, 1 - a * a);
//
//    J[2] = 1.0;
//    
//    J[3] = 0;
//
//}

double g_1(double u, double y) {
    double a = 5.0;
    double f = pow(10, 23);
    //double u = log(z / k);
    double dudt = -A * sqrt(B * exp(50.0 * u) + C * (1 - cos(y))) + D * sin(y) * exp(-50.0 * u);
    return dudt;
}
void Jacoby_Matrics_g_1(vector<double>& J, double u, double y) {
    double a = 5.0, f = pow(10, 23);
    J[0] = -A * B * 50.0 / 2.0 * exp(50.0 * u) / sqrt(B * exp(50 * u) + C * (1 - cos(y)))-50.0*D*sin(y)*exp(-50.0*u);

    J[1] = -A*C*sin(y)/(2*sqrt(B*exp(50*u)+C*(1-cos(y))))+D*cos(y)*exp(-50*u);

    J[2] = k*exp(u);

    J[3] = 0;

}

void Matrix_reverse(vector<double>& M) {
    if (M[0] * M[3] - M[1] * M[2] == 0) {
        cout << "no reverse Matrix";
    }
    else {
        double det_M = M[0] * M[3] - M[1] * M[2];
        vector<double> R = {
    M[3] / det_M,  -M[1] / det_M,
    -M[2] / det_M,  M[0] / det_M
        };
        M = R;
    }
}

void eq1(vector<double>& E1, vector<double>& J, double z, double y) {
    Jacoby_Matrics_g_1(J, z, y); //менять здеся
    E1[0] = 1 - h * g * J[0];
    E1[1] = 0 - h * g * J[1];
    E1[2] = 0 - h * g * J[2];
    E1[3] = 1-h*g*J[3];
}

//double integration(double z, double y, double I0) {
//    double a = 2*sqrt(24* 3.14159*G)*pow(10, 23);
//    double dr = pow(10, -3);
//    double f= sqrt(z * z + 2 - cos(y)) / (1 + L_4 * z * z);
//    return I0*pow(2.71828, -a * dr * f);
//}

double integration(double z, double y, double I0) {
    double a = 2*sqrt(24* 3.14159*G)*pow(10, 23);
    double dr = pow(10, -3);
    double f= sqrt(z * z + 2 - cos(y)) / (1 + L_4 * z * z);
    return I0*(1 + f * dr);
}
double kinetic_integration(double z, double y, double I0) {
    double a = 2 * sqrt(24 * 3.14159 * G) * pow(10, 23);
    double dr = pow(10, -3);
    double f = -sqrt(z * z + 2 - cos(y));
    return I0 * pow(2.71828, a * dr * f);
}
//void take_info( vector<double>& f, vector<double>& t, vector <double>& Y, double time) {
//    double buffer = a0 * sqrt(Y[0] * Y[0] + 2 - cos(Y[1])) / (1 + L_4 * Y[0] * Y[0]);
//    f.push_back(buffer);
//    t.push_back(time);
//}
//
//double integrate_adaptive_trapezoidal(const vector<double>& t, const vector<double>& f) {
//    double integral = 0.0;
//    for (size_t i = 1; i < t.size(); ++i) {
//        integral += 0.5 * (t[i] - t[i - 1]) * (f[i] + f[i - 1]);
//    }
//    return exp(integral);
