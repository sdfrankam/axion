// ConsoleApplication15.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

int main()
{
    double L_4 = pow(0.05, 4);
    ifstream csv_y, csv_z;
    ofstream rho_X, rho_V;
    csv_y.open("teta_q.txt");
    csv_z.open("dteta_q.txt");
    rho_X.open("rho_x.txt");
    rho_V.open("rho_v.txt");
    string line1, line2;
    int i = 0;
    double r0;
    while (getline(csv_z, line1) && getline(csv_y, line2)) {
        double z = stod(line1);
        double y = stod(line2);
        //КВИНТЕССЕНЦИЯ
        //r0 = z * z / 2 + pow(0.05, 4) * (1 - cos(y));
        //rho_X << ( z * z / (2.0*r0)) << endl;
        //rho_V << pow(0.05, 4) * (1 - cos(y))/r0 << endl;

        //ТАХИОНЫ
        //double L_4 = pow(0.05, 4);
        //r0 = L_4 * (1 - cos(y)) * (1+L_4*z*z/ (2*(1 - L_4 * z * z)));
        //rho_X << (L_4*L_4*z*z*(1-cos(y))/(2*(1-L_4*z*z)))/r0 << endl;
        //rho_V << L_4 *(1-cos(y))/ r0 << endl;

        //ФАНТОМНЫЕ ПОЛЯ
        //r0 = -z * z / 2 + pow(0.05, 4) * (1 - cos(y));
        //rho_X << ( -z * z / (2.0*r0)) << endl;
        //rho_V << pow(0.05, 4) * (1 - cos(y))/r0 << endl;

        //a>>1, b=0
        double a = 5.0;
        r0 = (2*a*a*a-a)*pow(z*z/2.0, a*a)+L_4*(1-cos(y));
        rho_X << (2 * a * a * a - a) * pow(z * z / 2.0, a * a)/r0 << endl;
        rho_V << L_4 *(1-cos(y))/ r0 << endl;


    }
    csv_y.close();
    csv_z.close();
    rho_X.close();
    rho_V.close();
}

