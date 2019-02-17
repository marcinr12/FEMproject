#pragma once
#include<vector>
#include<iostream>

using namespace std;

vector <vector <double>> multiplyVectorByTransparent(vector <double> v);
vector <double> multiplyVectorByNumber(vector <double> v, double d);
vector <vector <double>> multiplyMatrixByNumber(vector <vector <double>> v, double d);
vector <double> multiplyMatrixByVectorVector(vector<vector<double>> m, vector <double> u);
vector <vector <double>> divideMatrix(vector <vector <double>> v, double d);
vector <vector <double>> sumMatrix(vector <vector <double>> v, vector <vector <double>> u);
vector <double> sumVector(vector <double> v, vector <double> u);
void printMatrix(vector <vector <double>> v);
void printVector(vector <double> v);