#pragma once
#include <iostream>
#include <typeinfo>
#include <malloc.h>
using namespace std;

class SimpleFloatArray {
private:
    int num_elems;
    float* ptr_to_data;
    void copy(const SimpleFloatArray& a) {
        for (int i = 0; i < num_elems; i++) {
            ptr_to_data[i] = a.ptr_to_data[i];
        }
    };
public:
    SimpleFloatArray(int n) {
        num_elems = n;
        ptr_to_data = new float[n];
    }
    //SimpleFloatArray(int n, char x) {
    //  if (x == 's') {
    //    num_elems = n;
    //    ptr_to_data = new float[n];
    //    for (int i = 0; i < n; i++) {
    //      cin << ptr_to_data[i];
    //    }
    //  }
    //}
    SimpleFloatArray() {
        num_elems = 0;
        ptr_to_data = 0; 
    }
    SimpleFloatArray(const SimpleFloatArray& a) {
        num_elems = a.num_elems;
        ptr_to_data = new float[num_elems];
        copy(a); 
    }
    SimpleFloatArray(float* a, int size) {
        num_elems = size;
        ptr_to_data = new float[size];
        for (int i = 0; i < size; i++) {
            ptr_to_data[i] = a[i];
        }
        
    }
    ~SimpleFloatArray() {
        delete[] ptr_to_data;
    }; 

    void setSize(int n) {
        float* new_arr = new float[n];
        for (int i = 0; i < n; i++) {
            if (i < num_elems) {
                new_arr[i] = ptr_to_data[i];
            }
            else {
                new_arr[i] = 0;
            }
        }
        delete[] ptr_to_data;
        ptr_to_data = new_arr;
        num_elems = n;
    };

    float& operator[ ](int i) {
        return ptr_to_data[i];
    }; 

    int numElems() {
        int num_elems_pub = num_elems;
        return num_elems_pub;
    };
    SimpleFloatArray& operator=(const SimpleFloatArray& rhs) {
        if (ptr_to_data != rhs.ptr_to_data) {
            setSize(rhs.num_elems);
            copy(rhs);
        }
        return *this;
    }

    SimpleFloatArray& operator=(float rhs) {
        float* p = ptr_to_data + num_elems;
        while (p > ptr_to_data) *--p = rhs;
        return *this;

    };
};