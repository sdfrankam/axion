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
    ofstream /*rho_X, rho_V*/ data, full_energy;
    csv_y.open("teta_q.csv");
    csv_z.open("dteta_q.csv");
    //rho_X.open("rho_x.csv");
    //rho_V.open("rho_v.csv");
    data.open("data.csv");
    full_energy.open("full_energy.csv");
    string line1, line2, line2_y, line2_t;
    int i = 0;
    double r0;
    while (getline(csv_z, line1) && getline(csv_y, line2)) {
        double z = stod(line1);
	int j=0;
	for (char p: line2){
	if(static_cast<int>(p)!=static_cast<int>(',')){
	line2_y+=p;
	j++;
	}
	else break;
	}
	line2_t=line2.erase(0, j+1);
        //cout << line2_t << endl << line2_y << endl;
	double y = stod(line2_y);
	double t=stod(line2_t);
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
        data << t << "," << (2 * a * a * a - a) * pow(z * z / 2.0, a * a)/r0 << "," << L_4 *(1-cos(y))/ r0 <<  endl;
        full_energy << t << "," << r0 << endl;
    }
    csv_y.close();
    csv_z.close();
    full_energy.close();
    data.close();
    //rho_X.close();
    //rho_V.close();
}

