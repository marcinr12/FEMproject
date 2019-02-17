#pragma once
#include<vector>
#include"Grid.h"
#include"Element.h"

int il();

vector <vector<double>> calculateMatrixHdV(Element, double);
vector <vector<double>> calculateMatrixHdS(Element e, double convection, double heightGrid, double lenghGrid);
vector <vector<double>> calculateOneSideMatrixHdS(double ksi1, double ksi2, double eta1, double eta2, double lengh, double convection);
vector <vector <double>> calculateMatrixC(Element e, double c, double ro);
vector <double> calculateVectorP(Element e, double convection, double heightGrid, double lenghGrid, double ambientTemperature);
vector<double> calculateOneSideVectorP(double ksi1, double ksi2, double eta1, double eta2, double lengh, double convection, double ambientTemperature);

vector <vector <double>> calculateGlobalMatrixHdV(Grid g, double conductivity);
vector<vector<double>> calculateGlobalMatrixHdS(Grid g, double alfa, double heightGrid, double lenghGrid);
vector<vector<double>> calculateGlobalMatrixC(Grid g, double c, double ro);
vector <double> calculateGlobalVectorP(Grid g, double convection, double heightGrid, double lenghGrid, double ambientTemperature);