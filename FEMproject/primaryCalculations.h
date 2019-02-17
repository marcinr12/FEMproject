#pragma once
#include"Element.h"
#include"Node.h"
#include"Grid.h"

vector <double> N1calc();
vector <double> N2calc();
vector <double> N3calc();
vector <double> N4calc();

vector <double> calculateNsForPC(double ksi, double eta);

vector <double> dN_dKsi1Calc();
vector <double> dN_dKsi2Calc();
vector <double> dN_dKsi3Calc();
vector <double> dN_dKsi4Calc();

vector <double> dN_dEta1Calc();
vector <double> dN_dEta2Calc();
vector <double> dN_dEta3Calc();
vector <double> dN_dEta4Calc();

vector <double> J_1_1Cal(Element e, vector<double> dN_dKsi1, vector<double> dN_dKsi2, vector<double> dN_dKsi3, vector<double> dN_dKsi4);
vector <double> J_1_2Cal(Element e, vector<double> dN_dKsi1, vector<double> dN_dKsi2, vector<double> dN_dKsi3, vector<double> dN_dKsi4);
vector <double> J_2_1Cal(Element e, vector<double> dN_dEta1, vector<double> dN_dEta2, vector<double> dN_dEta3, vector<double> dN_dEta4);
vector <double> J_2_2Cal(Element e, vector<double> dN_dEta1, vector<double> dN_dEta2, vector<double> dN_dEta3, vector<double> dN_dEta4);

vector <double> detJCal(vector <double> J_1_1, vector <double> J_1_2, vector <double> J_2_1, vector <double> J_2_2);

vector < double> J_1_1_1Cal(vector<double>J_2_2, vector <double> detJ);
vector < double> J_1_1_2Cal(vector<double>J_1_2, vector <double> detJ);
vector < double> J_1_2_1Cal(vector<double>J_2_1, vector <double> detJ);
vector < double> J_1_2_2Cal(vector<double>J_1_1, vector <double> detJ);

vector <double> pc_dN_d(double j1, double j2, double dN_dK1, double dN_dK2, double dN_dK3, double dN_dK4, double dN_dE1, double dN_dE2, double dN_dE3, double dN_dE4);