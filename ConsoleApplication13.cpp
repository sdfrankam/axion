#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "fields.h"

//внесла изменения
void Step_1(vector<double>& E1, vector< double>& K1, double z, double y) {

    Matrix_reverse(E1);
    K1[0] = h * (E1[0] * g_1(z, y) + E1[1] * z);
    K1[1] = h * (E1[2] * g_1(z, y) + E1[3] * z);
}
void Step_2(vector< double>& K2, vector< double>& K1, vector<double>& E1, double z, double y) {
    vector<double> V1 = { (g_1(z + b*K1[0], y + b*K1[1])), (z +b* K1[0])};
    K2[0] =  h* (E1[0] * V1[0] + E1[1] * V1[1]);
    K2[1] =  h* (E1[2] * V1[0] + E1[3] * V1[1]);
}
void Step_3(vector<double>& Y_n, vector<double>& K1, vector<double>& K2, double z, double y) {

    Y_n[0] = z + p1 * K1[0] +p2 * K2[0];
    Y_n[1] = y + p1 * K1[1] +p2 * K2[1];
}



int main()
{
    ofstream csv_teta, csv_dteta, csv_time, csv_integral, csv_time_2;
    csv_teta.open("teta_q.csv");
    csv_dteta.open("dteta_q.csv");
    csv_time.open("time.txt");
    csv_integral.open("integral_kin.txt");
    csv_time_2.open("time_for_integral.txt");
    vector< double> Y_n = { -10.0, 3.1 };
    vector< double> K1 = { 0,0 }, K2 = { 0,0 }, E1{ 0,0,0,0 }, J = { 0,0,0,0 };
    vector<double> f, t;
    double I0 = 1;
    for (int i = 0; i < 3.5*pow(10, 6); i++) {
        double  z = Y_n[0], y = Y_n[1];
        eq1(E1, J, z, y);
        Step_1(E1, K1, z, y);
        Step_2(K2, K1, E1, z, y);
        Step_3(Y_n, K1, K2, z, y);
        //if (z>=400 || z<=-400) {
        //    cerr << "z is to big";
        //}
        csv_teta << Y_n[1] << "," << h*i << endl;
        csv_dteta << k*exp(Y_n[0]) << endl;
        if(h*i)
        //cout << Y_n[0] * Y_n[0] / 2 + Y_n[1] * Y_n[1] / 2 << endl;
        csv_time << h * i << endl;
        /*csv_integral << kinetic_integration(Y_n[0], Y_n[1], I0) << endl;
        I0 = integration(Y_n[0], Y_n[1], I0);
        csv_time_2 << 4 * i/3.124 << endl;*/
       //cout << "Step " << i << ": y = " << Y_n[1]
       //     << ", z = " << Y_n[0]
       //     << ", E = " << 0.5 * (Y_n[0] * Y_n[0] + Y_n[1] * Y_n[1])
       //     << endl;
    }
    cout << pow(3, 0) << endl;
    csv_teta.close();
    csv_time.close();
    csv_integral.close();
    csv_time_2.close();


}
