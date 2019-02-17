#pragma once

#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <string>

#include"Element.h"
#include"Node.h"
#include"GridSize.h"
#include"Grid.h"
#include"primaryCalculations.h"
#include"matrixOperations.h"
#include"MESmatrices.h"
#include"jacobiMethod.h"
#include "Saver.h"

using namespace std;

void tc1();
void tc2();
void printCSV(vector<double> vectorP, int nH, int nL, string s, string iteracja)
{
	std::ofstream plikW;
	plikW.open("Element" + s + ".csv", ios::out | ios::app);// , std::ios::in | std::ios::out);
	plikW << endl << endl;
	plikW << "Iteracja: " << iteracja << endl;
	for (int i = 0; i < nH; i++)
	{
		for (int j = 0; j < nL; j++)
		{
			plikW << vectorP[i * nL + j] << ';';
			//cout << vectorP[i * nL +  j] << '\t';
		}
		plikW << endl;
		//cout << endl;
	}
}

vector<double> gauss(vector< vector<double> > matrix, vector<double> coefficients) 
{

	int n = matrix.size();

	vector<double> line(n + 1, 0);
	vector< vector<double> > A(n, line);

	// Read input data
	for (int i = 0; i<n; i++) {
		for (int j = 0; j<n; j++) {
			A[i][j] = matrix[i][j];
		}
	}

	for (int i = 0; i<n; i++) {
		A[i][n] = coefficients[i];
	}

	for (int i = 0; i<n; i++) {
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = abs(A[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k<n + 1; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k<n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j<n + 1; j++) {
				if (i == j) {
					A[k][j] = 0;
				}
				else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = A[i][n] / A[i][i];
		for (int k = i - 1; k >= 0; k--) {
			A[k][n] -= A[k][i] * x[i];
		}
	}
	return x;
}

double fun(double k, double e)
{
	return (5+3*k)*(1+e) + (1+e);
}

int main()
{

	//cout.precision(5);

	//GridSize *gs = new GridSize();
	//getData(gs);						//odczyt z pliku
	//gs->getData();

	//cout << "H: " << gs->H << endl << "L: " << gs->L << endl << "nH: " << gs->nH << endl << "nL: " << gs->nL << endl;
	//cout << endl;


	//Grid *grid = new Grid;
								//40
								//30
								//6
								//4
	//grid->createGrid(gs->nH, gs->nL, gs->H, gs->L, 100);

	//grid->printGrid();		//wypisanie wez³ów zgodnie z numeracj¹

	//Node(x, y, coefficient)
	//Node *n1 = new Node(0, 0.0333, 1);
	//Node *n2 = new Node(0.0333, 0.0333, 1);
	//Node *n3 = new Node(0.0333, 0.0666, 1);
	//Node *n4 = new Node(0, 0.0666, 1);

	//Element *e = new Element(n1, n2, n3, n4);
	//cout << e->nodes[0]->x << endl;
	//cout << n1->x << endl;

	//Element *e = new Element(0, 0, 0.025, 0, 0.025, 0.025, 0, 0.025, 1, 1, 1, 1);
	//cout << e->nodes[1]->x << endl;

	//cout << "MATRIX H" << endl;
	//vector <vector <double>> matrixH = calculateMatrixHdV(*e, 25);
	//printMatrix2d(matrixH);

	//cout << "MATRIX C" << endl;
	//vector <vector <double>> matrixC = calculateMatrixC(*e, 700, 7800);
	//printMatrix2d(matrixC);

	//cout << "VECTOR P" << endl;
	//printMatrix1d(calculateVectorP(*e));

	
	//////////////////////////////////////
	//-----initial data-----
	//////////////////////////////////////
	//double dt = 50;
	//double ambientTemperature = 1200;
	//double conductivity = 25;
	//double convection = 300;
	//double c = 700;
	//double ro = 7800;
	//vector<double> t0(grid->nodes.size());
	//for (int i = 0; i < t0.size(); i++)
	//	t0[i] = 100;
	//////////////////////////////////////

	////cout << "GLOBAL MATRIX HdV" << endl;
	//vector <vector <double>> hgv = calculateGlobalMatrixHdV(*grid, conductivity);
	////printMatrix2d(hgv);

	/*cout << "GLOBAL MATRIX HdS" << endl;
	vector <vector <double>> hgs = calculateGlobalMatrixHdS(*grid, convection, gs->H, gs->L);
	printMatrix2d(hgs);*/

	////cout << "GLOBAL MATRIX HdV + HdS" << endl;
	//vector <vector <double>> hgvs = sumMatrix2d(hgs, hgv);
	////printMatrix2d(hgvs);

	/*cout << "GLOBAL MATRIX C" << endl;
	vector <vector <double>> cg = calculateGlobalMatrixC(*grid, c, ro);
	printMatrix(cg);*/

	//cout << "[H] = [H] + [C]/dT" << endl;
	//vector <vector <double>> hcdt = sumMatrix2d(hgvs, devideMatrix2d(cg, dt));
	//printMatrix2d(hcdt);

	////cout << "[C]/dT" << endl;
	//vector <vector <double>> cdt = devideMatrix2d(cg, dt);
	//////printMatrix2d(cdt);

	//cout << "GLOBAL VECTOR P" << endl;
	//vector <double> pg = calculateGlobalVectorP(*grid, convection, gs->H, gs->L, ambientTemperature);
	//printMatrix1d(pg);

	////cout << "{[C]/dT} * {T0}" << endl;
	//vector <double> cdtt0 = multiplyMatrixByVectorVector(cdt, t0);
	////printMatrix1d(cdtt0);

	//cout << "([{P} + {[C]/dT}*{T0})" << endl;
	//vector <double> Pcdtt0 = sumMatrix1d(pg, cdtt0);
	//printMatrix1d(Pcdtt0);

	//cout << "{T0}" << endl;
	//printMatrix1d(t0);

	//for (int i = 0; i < Pcdtt0.size(); i++)
	//	Pcdtt0[i] *= -1.0;

	////			wsp => H	roz => {P}		rzad => int ile iteracji (dok³adnosc rozwiazan)
	//cout << "METODA JACOBIEGO" << endl;
	//printMatrix1d(metodaJacobiego(hcdt, Pcdtt0, t0, 10));

	//tc1();
	tc2();
	

	//cout << endl;
	

	system("PAUSE");
	return 0;
}

void tc1()
{
	GridSize *gs = new GridSize();
	//getData(gs);						//odczyt z pliku
	gs->getData();
	gs->H = 0.1;
	gs->L = 0.1;
	gs->nH = 4;
	gs->nL = 4;

	cout << "H: " << gs->H << endl << "L: " << gs->L << endl << "nH: " << gs->nH << endl << "nL: " << gs->nL << endl;
	cout << endl;

	Grid *grid = new Grid;

	// = createGrid(*g, gs->nH, gs->nL, gs->H, gs->L, 100);	//tworzy grid (tworzy i numeruje odpowiednio wez³y)
	grid->createGrid(gs->nH, gs->nL, gs->H, gs->L);
	//###############################################################################
	//-----------------------TEST CASE 1---------------------------------------------
	//###############################################################################

	//////////////////////////////////////
	//-----initial data-----
	//////////////////////////////////////
	double dt = 50;
	double ambientTemperature = 1200;
	double conductivity = 25;
	double convection = 300;
	double c = 700;
	double ro = 7800;
	vector<double> t0(grid->nodes.size());
	for (int i = 0; i < t0.size(); i++)
		t0[i] = 100;
	//////////////////////////////////////
	double maxH = 0;
	double maxL = 0;
	for (int i = 0; i < grid->elements.size(); i++)
	{
		for (int j = 0; j < grid->elements[i].nodes.size(); j++)
		{
			if (grid->elements[i].nodes[j]->y > maxH)
				maxH = grid->elements[i].nodes[j]->y;
			if (grid->elements[i].nodes[j]->x > maxL)
				maxL = grid->elements[i].nodes[j]->x;
		}
	}
	std::ofstream plikM;
	plikM.open("MinMaxTC1.csv", ios::out | ios::app);// , std::ios::in | std::ios::out);
	clock_t startM = clock();
	vector <vector <double>> hgv = calculateGlobalMatrixHdV(*grid, conductivity);
	vector <vector <double>> hgs = calculateGlobalMatrixHdS(*grid, convection, maxH, maxL);
	vector <vector <double>> hgvs = sumMatrix(hgs, hgv);
	vector <vector <double>> cg = calculateGlobalMatrixC(*grid, c, ro);
	vector <vector <double>> hcdt = sumMatrix(hgvs, divideMatrix(cg, dt));
	vector <vector <double>> cdt = divideMatrix(cg, dt);
	std::ofstream plikW;
	plikW.open("WynikiTC1.txt");// , std::ios::in | std::ios::out);
	cout << "Time matrices: " << (clock() - startM) / 1000.0 << "[s]" << endl;
	plikW << "Time matrices: " << (clock() - startM) / 1000.0 << "[s]" << endl;

	cout << "ELIMINACJA GAUSSA" << endl;
	plikW << "ELIMINACJA GAUSSA" << endl;
	int iteracje = 10;
	for (int i = 0; i < iteracje; i++)
	{
		clock_t start = clock();
		vector <double> pg = calculateGlobalVectorP(*grid, convection, maxH, maxL, ambientTemperature);
		vector <double> cdtt0 = multiplyMatrixByVectorVector(cdt, t0);
		vector <double> Pcdtt0 = sumVector(pg, cdtt0);


		//for (int j = 0; j < Pcdtt0.size(); j++)
		//Pcdtt0[j] *= -1.0;

		//			wsp => H	roz => {P}		rzad => int ile iteracji (dok³adnosc rozwiazan)
		clock_t startJ = clock();
		vector<double> wyniki = gauss(hcdt, Pcdtt0);
		//vector<double> wyniki = metodaJacobiego(hcdt, Pcdtt0, t0, 7);
		cout << "Iteracja: " << i + 1 << '\t' << "Gauss: " << (clock() - startJ) / 1000.0 << "[s]" << '\t';
		plikW << "Iteracja: " << i + 1 << '\t' << "Gauss: " << (clock() - startJ) / 1000.0 << "[s]" << '\t';
		//printVector(wyniki);
		double max = *(std::max_element(wyniki.begin(), wyniki.end()));
		double min = *(std::min_element(wyniki.begin(), wyniki.end()));
		cout << "Max value: " << max << '\t' << "Min value: " << min << '\t' << "Total time: " << (clock() - start) / 1000.0 << "[s]" << endl;
		plikW << "Max value: " << max << '\t' << "Min value: " << min << '\t' << "Total time: " << (clock() - start) / 1000.0 << "[s]" << endl;
		printCSV(wyniki, gs->nH, gs->nL, "TC" + to_string(1), to_string(i + 1));
		plikM << min << ';' << max << endl;
		//printVector(wyniki);
		/*cout << "SIZE" << t0.size() << endl;*/

		for (int j = 0; j < t0.size(); j++)
			t0[j] = wyniki[j];
	}
	plikM.close();
	plikW.close();
}

void tc2()
{

	GridSize *gs = new GridSize();
	//getData(gs);						//odczyt z pliku
	gs->getData();
	gs->H = 0.1;
	gs->L = 0.1;
	gs->nH = 31;
	gs->nL = 31;

	cout << "H: " << gs->H << endl << "L: " << gs->L << endl << "nH: " << gs->nH << endl << "nL: " << gs->nL << endl;
	cout << endl;

	Grid *grid = new Grid;

	// = createGrid(*g, gs->nH, gs->nL, gs->H, gs->L, 100);	//tworzy grid (tworzy i numeruje odpowiednio wez³y)
	grid->createGrid(gs->nH, gs->nL, gs->H, gs->L);
	//###############################################################################
	//-----------------------TEST CASE 2---------------------------------------------
	//###############################################################################

	//////////////////////////////////////
	//-----initial data-----
	//////////////////////////////////////
	double dt = 1;
	double ambientTemperature = 1200;
	double conductivity = 25;
	double convection = 300;
	double c = 700;
	double ro = 7800;
	vector<double> t0(grid->nodes.size());
	for (int i = 0; i < t0.size(); i++)
		t0[i] = 100;
	//////////////////////////////////////
	double maxH = 0;
	double maxL = 0;
	for (int i = 0; i < grid->elements.size(); i++)
	{
		for (int j = 0; j < grid->elements[i].nodes.size(); j++)
		{
			if (grid->elements[i].nodes[j]->y > maxH)
				maxH = grid->elements[i].nodes[j]->y;
			if (grid->elements[i].nodes[j]->x > maxL)
				maxL = grid->elements[i].nodes[j]->x;
		}
	}
	std::ofstream plikM;
	plikM.open("MinMaxTC2.csv", ios::out | ios::app);// , std::ios::in | std::ios::out);
	clock_t startM = clock();
	vector <vector <double>> hgv = calculateGlobalMatrixHdV(*grid, conductivity);
	vector <vector <double>> hgs = calculateGlobalMatrixHdS(*grid, convection, maxH, maxL);
	vector <vector <double>> hgvs = sumMatrix(hgs, hgv);
	vector <vector <double>> cg = calculateGlobalMatrixC(*grid, c, ro);
	vector <vector <double>> hcdt = sumMatrix(hgvs, divideMatrix(cg, dt));
	vector <vector <double>> cdt = divideMatrix(cg, dt);
	std::ofstream plikW;
	plikW.open("WynikiTC2.txt");// , std::ios::in | std::ios::out);
	cout << "Time matrices: " << (clock() - startM) / 1000.0 << "[s]" << endl;
	plikW << "Time matrices: " << (clock() - startM) / 1000.0 << "[s]" << endl;

	cout << "ELIMINACJA GAUSSA" << endl;
	plikW << "ELIMINACJA GAUSSA" << endl;
	int iteracje = 100;
	Saver::save_to_file("Sim.txt", iteracje);
	for (int i = 0; i < iteracje; i++)
	{
		clock_t start = clock();
		vector <double> pg = calculateGlobalVectorP(*grid, convection, maxH, maxL, ambientTemperature);
		vector <double> cdtt0 = multiplyMatrixByVectorVector(cdt, t0);
		vector <double> Pcdtt0 = sumVector(pg, cdtt0);


		//for (int j = 0; j < Pcdtt0.size(); j++)
		//Pcdtt0[j] *= -1.0;

		//			wsp => H	roz => {P}		rzad => int ile iteracji (dok³adnosc rozwiazan)
		clock_t startJ = clock();
		vector<double> wyniki = gauss(hcdt, Pcdtt0);
		//vector<double> wyniki = metodaJacobiego(hcdt, Pcdtt0, t0, 7);
		cout << "Iteracja: " << i + 1 << '\t' << "Gauss: " << (clock() - startJ) / 1000.0 << "[s]" << '\t';
		plikW << "Iteracja: " << i + 1 << '\t' << "Gauss: " << (clock() - startJ) / 1000.0 << "[s]" << '\t';
		//printVector(wyniki);
		double max = *(std::max_element(wyniki.begin(), wyniki.end()));
		double min = *(std::min_element(wyniki.begin(), wyniki.end()));
		cout << "Max value: " << max << '\t' << "Min value: " << min << '\t' << "Total time: " << (clock() - start) / 1000.0 << "[s]" << endl;
		plikW << "Max value: " << max << '\t' << "Min value: " << min << '\t' << "Total time: " << (clock() - start) / 1000.0 << "[s]" << endl;
		printCSV(wyniki, gs->nH, gs->nL, "TC" + to_string(1), to_string(i + 1));
		plikM << min << ';' << max << endl;
		//printVector(wyniki);
		/*cout << "SIZE" << t0.size() << endl;*/
		Saver::save_to_file("Sim.txt", wyniki);
		for (int j = 0; j < t0.size(); j++)
			t0[j] = wyniki[j];
	}
	plikM.close();
	plikW.close();
}
