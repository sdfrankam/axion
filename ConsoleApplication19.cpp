#include <iostream>
#include <typeinfo>
#include <malloc.h>
#include "Header.h"
using namespace std;

int main()
{
    int n;
    cout << "Print size of array : " << endl;
    cin >> n;
    float* buffer_array = new float[n];
    for (int i = 0; i < n; i++) {
        float buffer_el;
        cout << "Print element of array n. " << i + 1 << endl;
        cin >> buffer_el;
        buffer_array[i] = buffer_el;

    }
    
    SimpleFloatArray new_arr(buffer_array, n);
    new_arr.setSize(7);
    cout << "size of array increase to 7 : " << endl;
    for (int i = 0; i < new_arr.numElems(); i++) {
        cout << new_arr[i] << " ";
    };
    cout << endl;
    new_arr.setSize(3);
    cout << "size of array decrease to 3 : " << endl;
    for (int i = 0; i < new_arr.numElems(); i++) {
        cout << new_arr[i] << " ";
    };
    cout << endl;
    new_arr.setSize(7);
    cout << "size of array increase to 7 again : " << endl;
    for (int i = 0; i < new_arr.numElems(); i++) {
        cout << new_arr[i] << " ";
    };
    cout << endl;
    cout << "numElems function " << new_arr.numElems() << endl;

    //new_arr.~SimpleFloatArray();

}