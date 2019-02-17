#pragma once

#include <vector>
using namespace std;

vector<vector<double>> podzielL(vector<vector<double>> tab);

vector<vector<double>> podzielD(vector<vector<double>> tab);

vector<vector<double>> podzielU(vector<vector<double>> tab);

vector<vector<double>> obliczN(vector<vector<double>> tab);

vector<vector<double>> obliczSume(vector<vector<double>> tab1, vector<vector<double>> tab2);

vector<vector<double>> obliczIloczyn(vector<vector<double>> tab1, vector<vector<double>> tab2);

vector<vector<double>> obliczPrzeciwna(vector<vector<double>> tab);

vector<double> metodaJacobiego(vector<vector<double>> wsp, vector<double> roz, vector <double> przyblizoneRozw, int rzad);

