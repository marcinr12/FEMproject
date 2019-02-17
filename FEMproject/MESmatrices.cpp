#include"MESmatrices.h"
#include"primaryCalculations.h"
#include"matrixOperations.h"

vector <vector <double>> calculateMatrixHdV(Element e, double k)
{
	double ksi1 = -1.0 / sqrt(3), eta1 = ksi1;
	double ksi2 = 1.0 / sqrt(3), eta2 = eta1;
	double ksi3 = ksi2, eta3 = ksi3;
	double ksi4 = -ksi2, eta4 = eta3;;

	vector <double>  dN_dKsi1 = dN_dKsi1Calc();
	vector <double>  dN_dKsi2 = dN_dKsi2Calc();
	vector <double>  dN_dKsi3 = dN_dKsi3Calc();
	vector <double>  dN_dKsi4 = dN_dKsi4Calc();

	vector <double>  dN_dEta1 = dN_dEta1Calc();
	vector <double>  dN_dEta2 = dN_dEta2Calc();
	vector <double>  dN_dEta3 = dN_dEta3Calc();
	vector <double>  dN_dEta4 = dN_dEta4Calc();

	vector <double> J_1_1 = J_1_1Cal(e, dN_dKsi1, dN_dKsi2, dN_dKsi3, dN_dKsi4);
	vector <double> J_1_2 = J_1_2Cal(e, dN_dKsi1, dN_dKsi2, dN_dKsi3, dN_dKsi4);
	vector <double> J_2_1 = J_2_1Cal(e, dN_dEta1, dN_dEta2, dN_dEta3, dN_dEta4);
	vector <double> J_2_2 = J_2_2Cal(e, dN_dEta1, dN_dEta2, dN_dEta3, dN_dEta4);

	vector <double> detJ = detJCal(J_1_1, J_1_2, J_2_1, J_2_2);

	vector <double> J_1_1_1 = J_1_1_1Cal(J_2_2, detJ);
	vector <double> J_1_1_2 = J_1_1_2Cal(J_1_2, detJ);
	vector <double> J_1_2_1 = J_1_2_1Cal(J_2_1, detJ);
	vector <double> J_1_2_2 = J_1_2_2Cal(J_1_1, detJ);

	vector <double> pc1xTabV = pc_dN_d(J_1_1_1[0], J_1_1_2[0], dN_dKsi1[0], dN_dKsi2[0], dN_dKsi3[0], dN_dKsi4[0], dN_dEta1[0], dN_dEta2[0], dN_dEta3[0], dN_dEta4[0]);
	vector <double> pc2xTabV = pc_dN_d(J_1_1_1[1], J_1_1_2[1], dN_dKsi1[1], dN_dKsi2[1], dN_dKsi3[1], dN_dKsi4[1], dN_dEta1[1], dN_dEta2[1], dN_dEta3[1], dN_dEta4[1]);
	vector <double> pc3xTabV = pc_dN_d(J_1_1_1[2], J_1_1_2[2], dN_dKsi1[2], dN_dKsi2[2], dN_dKsi3[2], dN_dKsi4[2], dN_dEta1[2], dN_dEta2[2], dN_dEta3[2], dN_dEta4[2]);
	vector <double> pc4xTabV = pc_dN_d(J_1_1_1[3], J_1_1_2[3], dN_dKsi1[3], dN_dKsi2[3], dN_dKsi3[3], dN_dKsi4[3], dN_dEta1[3], dN_dEta2[3], dN_dEta3[3], dN_dEta4[3]);

	vector <double> pc1yTabV = pc_dN_d(J_1_2_1[0], J_1_2_2[0], dN_dKsi1[0], dN_dKsi2[0], dN_dKsi3[0], dN_dKsi4[0], dN_dEta1[0], dN_dEta2[0], dN_dEta3[0], dN_dEta4[0]);
	vector <double> pc2yTabV = pc_dN_d(J_1_2_1[1], J_1_2_2[1], dN_dKsi1[1], dN_dKsi2[1], dN_dKsi3[1], dN_dKsi4[1], dN_dEta1[1], dN_dEta2[1], dN_dEta3[1], dN_dEta4[1]);
	vector <double> pc3yTabV = pc_dN_d(J_1_2_1[2], J_1_2_2[2], dN_dKsi1[2], dN_dKsi2[2], dN_dKsi3[2], dN_dKsi4[2], dN_dEta1[2], dN_dEta2[2], dN_dEta3[2], dN_dEta4[2]);
	vector <double> pc4yTabV = pc_dN_d(J_1_2_1[3], J_1_2_2[3], dN_dKsi1[3], dN_dKsi2[3], dN_dKsi3[3], dN_dKsi4[3], dN_dEta1[3], dN_dEta2[3], dN_dEta3[3], dN_dEta4[3]);

	vector <vector<double>> pc1xTV = multiplyVectorByTransparent(pc1xTabV);
	vector <vector<double>> pc2xTV = multiplyVectorByTransparent(pc2xTabV);
	vector <vector<double>> pc3xTV = multiplyVectorByTransparent(pc3xTabV);
	vector <vector<double>> pc4xTV = multiplyVectorByTransparent(pc4xTabV);

	vector <vector<double>> pc1yTV = multiplyVectorByTransparent(pc1yTabV);
	vector <vector<double>> pc2yTV = multiplyVectorByTransparent(pc2yTabV);
	vector <vector<double>> pc3yTV = multiplyVectorByTransparent(pc3yTabV);
	vector <vector<double>> pc4yTV = multiplyVectorByTransparent(pc4yTabV);

	vector <vector<double>> pc1xTDetV = multiplyMatrixByNumber(pc1xTV, detJ[0]);
	vector <vector<double>> pc2xTDetV = multiplyMatrixByNumber(pc2xTV, detJ[1]);
	vector <vector<double>> pc3xTDetV = multiplyMatrixByNumber(pc3xTV, detJ[2]);
	vector <vector<double>> pc4xTDetV = multiplyMatrixByNumber(pc4xTV, detJ[3]);

	vector <vector<double>> pc1yTDetV = multiplyMatrixByNumber(pc1yTV, detJ[0]);
	vector <vector<double>> pc2yTDetV = multiplyMatrixByNumber(pc2yTV, detJ[1]);
	vector <vector<double>> pc3yTDetV = multiplyMatrixByNumber(pc3yTV, detJ[2]);
	vector <vector<double>> pc4yTDetV = multiplyMatrixByNumber(pc4yTV, detJ[3]);

	vector <vector<double>> sum1V = sumMatrix(pc1xTDetV, pc1yTDetV);
	vector <vector<double>> sum2V = sumMatrix(pc2xTDetV, pc2yTDetV);
	vector <vector<double>> sum3V = sumMatrix(pc3xTDetV, pc3yTDetV);
	vector <vector<double>> sum4V = sumMatrix(pc4xTDetV, pc4yTDetV);

	vector <vector<double>> sum1kV = multiplyMatrixByNumber(sum1V, k);
	vector <vector<double>> sum2kV = multiplyMatrixByNumber(sum2V, k);
	vector <vector<double>> sum3kV = multiplyMatrixByNumber(sum3V, k);
	vector <vector<double>> sum4kV = multiplyMatrixByNumber(sum4V, k);

	vector <vector<double>> matrixH = sumMatrix(sumMatrix(sum1kV, sum2kV), sumMatrix(sum3kV, sum4kV));

	return matrixH;
}

vector <vector<double>> calculateMatrixHdS(Element e, double convection, double heightGrid, double lengthGrid)
{
	double lenghSite1 = sqrt(pow((e.nodes[1]->x - e.nodes[0]->x), 2) + pow((e.nodes[1]->y - e.nodes[0]->y), 2));		//beetwen nodes 0 - 1
	double lenghSite2 = sqrt(pow((e.nodes[2]->x - e.nodes[1]->x), 2) + pow((e.nodes[2]->y - e.nodes[1]->y), 2));		//beetwen nodes 1 - 2
	double lenghSite3 = sqrt(pow((e.nodes[3]->x - e.nodes[2]->x), 2) + pow((e.nodes[3]->y - e.nodes[2]->y), 2));		//beetwen nodes 3 - 2
	double lenghSite4 = sqrt(pow((e.nodes[0]->x - e.nodes[1]->x), 2) + pow((e.nodes[0]->y - e.nodes[1]->y), 2));		//beetwen nodes 1 - 0

	bool onSide1 = false;
	bool onSide2 = false;
	bool onSide3 = false;
	bool onSide4 = false;
	

	if (e.nodes[0]->y == 0 && e.nodes[1]->y == 0) 
		onSide1 = true;
	
	if (e.nodes[1]->x == lengthGrid && e.nodes[2]->x == lengthGrid)
		onSide2 = true;
	
		
	if (e.nodes[2]->y == heightGrid && e.nodes[3]->y == heightGrid)
		onSide3 = true;
	
		
	if (e.nodes[3]->x == 0 && e.nodes[0]->x == 0)
		onSide4 = true;
	
	vector<vector<double>> sum1(4);
	vector<vector<double>> sum2(4);
	vector<vector<double>> sum3(4);
	vector<vector<double>> sum4(4);
	for (int i = 0; i < sum1.size(); i++)
	{
		sum1[i].resize(4);
		sum2[i].resize(4);
		sum3[i].resize(4);
		sum4[i].resize(4);
	}
	for (int i = 0; i < sum1.size(); i++)
		for (int j = 0; j < sum1[0].size(); j++)
		{
			sum1[i][j] = 0;
			sum2[i][j] = 0;
			sum3[i][j] = 0;
			sum4[i][j] = 0;
		}

	if (onSide1)
	{
		double ksi1 = -1 / (sqrt(3));
		double ksi2 = 1 / (sqrt(3));
		double eta1 = -1;
		double eta2 = -1;

		sum1 = calculateOneSideMatrixHdS(ksi1, ksi2, eta1, eta2, lenghSite1, convection);
	}

	if (onSide2)
	{
		double ksi1 = 1;
		double ksi2 = 1;
		double eta1 = -1 / (sqrt(3));
		double eta2 = 1 / (sqrt(3));

		sum2 = calculateOneSideMatrixHdS(ksi1, ksi2, eta1, eta2, lenghSite2, convection);
	}

	if (onSide3)
	{
		double ksi1 = 1 / (sqrt(3));
		double ksi2 = -1 / (sqrt(3));
		double eta1 = 1;
		double eta2 = 1;

		sum3 = calculateOneSideMatrixHdS(ksi1, ksi2, eta1, eta2, lenghSite3, convection);
	}

	if (onSide4)
	{
		double ksi1 = -1;
		double ksi2 = -1;
		double eta1 = 1 / (sqrt(3));
		double eta2 = -1 / (sqrt(3));

		sum4 = calculateOneSideMatrixHdS(ksi1, ksi2, eta1, eta2, lenghSite4, convection);
	}
	vector<vector<double>> matrixHdS = sumMatrix(sumMatrix(sum1, sum2), sumMatrix(sum3, sum4));

	return matrixHdS;
}

vector <vector<double>> calculateOneSideMatrixHdS(double ksi1, double ksi2, double eta1, double eta2, double length, double convection)
{
	double detJ = length / 2.0;

	vector <double> pc1v = calculateNsForPC(ksi1, eta1);
	vector <double> pc2v = calculateNsForPC(ksi2, eta2);

	vector <vector <double>> pc1m = multiplyVectorByTransparent(pc1v);
	vector <vector <double>> pc2m = multiplyVectorByTransparent(pc2v);

	vector <vector <double>> pc1ma = multiplyMatrixByNumber(pc1m, convection);
	vector <vector <double>> pc2ma = multiplyMatrixByNumber(pc2m, convection);

	vector <vector <double>> sum = sumMatrix(pc1ma, pc2ma);
	vector <vector <double>> sumDetJ = multiplyMatrixByNumber(sum, detJ);

	return sumDetJ;
}

vector <vector <double>> calculateMatrixC(Element e, double c, double ro)
{
	vector <double>  dN_dKsi1 = dN_dKsi1Calc();
	vector <double>  dN_dKsi2 = dN_dKsi2Calc();
	vector <double>  dN_dKsi3 = dN_dKsi3Calc();
	vector <double>  dN_dKsi4 = dN_dKsi4Calc();

	vector <double>  dN_dEta1 = dN_dEta1Calc();
	vector <double>  dN_dEta2 = dN_dEta2Calc();
	vector <double>  dN_dEta3 = dN_dEta3Calc();
	vector <double>  dN_dEta4 = dN_dEta4Calc();

	vector <double> J_1_1 = J_1_1Cal(e, dN_dKsi1, dN_dKsi2, dN_dKsi3, dN_dKsi4);
	vector <double> J_1_2 = J_1_2Cal(e, dN_dKsi1, dN_dKsi2, dN_dKsi3, dN_dKsi4);
	vector <double> J_2_1 = J_2_1Cal(e, dN_dEta1, dN_dEta2, dN_dEta3, dN_dEta4);
	vector <double> J_2_2 = J_2_2Cal(e, dN_dEta1, dN_dEta2, dN_dEta3, dN_dEta4);

	vector <double> detJ = detJCal(J_1_1, J_1_2, J_2_1, J_2_2);

	vector <double> N1 = N1calc();
	vector <double> N2 = N2calc();
	vector <double> N3 = N3calc();
	vector <double> N4 = N4calc();

	vector<vector<double>> pc1(4);
	vector<vector<double>> pc2(4);
	vector<vector<double>> pc3(4);
	vector<vector<double>> pc4(4);
	for (int i = 0; i < pc1.size(); i++)
	{
		pc1[i].resize(4);
		pc2[i].resize(4);
		pc3[i].resize(4);
		pc4[i].resize(4);
	}

	for (int i = 0; i < pc1.size(); i++)
	{
		for (int j = 0; j < pc1[0].size(); j++)
		{
			pc1[i][j] = N1[i] * N1[j] * detJ[j] * c * ro;
			pc2[i][j] = N2[i] * N2[j] * detJ[j] * c * ro;
			pc3[i][j] = N3[i] * N3[j] * detJ[j] * c * ro;
			pc4[i][j] = N4[i] * N4[j] * detJ[j] * c * ro;

		}
	}

	vector <vector <double>> matrixC = sumMatrix(sumMatrix(pc1, pc2), sumMatrix(pc3, pc4));

	return matrixC;
}

vector <double> calculateVectorP(Element e, double convection, double heightGrid, double lenghGrid, double ambientTemperature)
{
	double lenghSite1 = sqrt(pow((e.nodes[1]->x - e.nodes[0]->x), 2) + pow((e.nodes[1]->y - e.nodes[0]->y), 2));		//beetwen nodes 0 - 1
	double lenghSite2 = sqrt(pow((e.nodes[2]->x - e.nodes[1]->x), 2) + pow((e.nodes[2]->y - e.nodes[1]->y), 2));		//beetwen nodes 1 - 2
	double lenghSite3 = sqrt(pow((e.nodes[3]->x - e.nodes[2]->x), 2) + pow((e.nodes[3]->y - e.nodes[2]->y), 2));		//beetwen nodes 3 - 2
	double lenghSite4 = sqrt(pow((e.nodes[0]->x - e.nodes[3]->x), 2) + pow((e.nodes[0]->y - e.nodes[3]->y), 2));		//beetwen nodes 3 - 0

	bool onSide1 = false;
	bool onSide2 = false;
	bool onSide3 = false;
	bool onSide4 = false;

	if (e.nodes[0]->y == 0 && e.nodes[1]->y == 0)
		onSide1 = true;
	if (e.nodes[1]->x == lenghGrid && e.nodes[2]->x == lenghGrid)
		onSide2 = true;
	if (e.nodes[2]->y == heightGrid && e.nodes[3]->y == heightGrid)
		onSide3 = true;
	if (e.nodes[3]->x == 0 && e.nodes[0]->x == 0)
		onSide4 = true;

	vector<double> sum1(4);
	vector<double> sum2(4);
	vector<double> sum3(4);
	vector<double> sum4(4);

	for (int i = 0; i < sum1.size(); i++)
	{
		sum1[i] = 0;
		sum2[i] = 0;
		sum3[i] = 0;
		sum4[i] = 0;
	}

	if (onSide1)
	{
		double ksi1 = -1 / (sqrt(3));
		double ksi2 = 1 / (sqrt(3));
		double eta1 = -1;
		double eta2 = -1;

		sum1 = calculateOneSideVectorP(ksi1, ksi2, eta1, eta2, lenghSite1, convection, ambientTemperature);
	}

	if (onSide2)
	{
		double ksi1 = 1;
		double ksi2 = 1;
		double eta1 = -1 / (sqrt(3));
		double eta2 = 1 / (sqrt(3));

		sum2 = calculateOneSideVectorP(ksi1, ksi2, eta1, eta2, lenghSite2, convection, ambientTemperature);
	}

	if (onSide3)
	{
		double ksi1 = 1 / (sqrt(3));
		double ksi2 = -1 / (sqrt(3));
		double eta1 = 1;
		double eta2 = 1;

		sum3 = calculateOneSideVectorP(ksi1, ksi2, eta1, eta2, lenghSite3, convection, ambientTemperature);
	}

	if (onSide4)
	{
		double ksi1 = -1;
		double ksi2 = -1;
		double eta1 = 1 / (sqrt(3));
		double eta2 = -1 / (sqrt(3));

		sum4 = calculateOneSideVectorP(ksi1, ksi2, eta1, eta2, lenghSite4, convection, ambientTemperature);
	}
	vector <double> vectorP = sumVector(sumVector(sum1, sum2), sumVector(sum3, sum4));

	return vectorP;

}

vector<double> calculateOneSideVectorP(double ksi1, double ksi2, double eta1, double eta2, double lengh, double convection, double ambientTemperature)
{
	double detJ = lengh / 2;

	vector <double> pc1v = calculateNsForPC(ksi1, eta1);
	vector <double> pc2v = calculateNsForPC(ksi2, eta2);

	vector <double> pc1vta = multiplyVectorByNumber(pc1v, convection * ambientTemperature);
	vector <double> pc2vta = multiplyVectorByNumber(pc2v, convection * ambientTemperature);

	vector <double> sum = sumVector(pc1vta, pc2vta);
	vector <double> sumtaDetJ = multiplyVectorByNumber(sum, detJ);

	return sumtaDetJ;
}

vector<vector<double>> calculateGlobalMatrixHdV(Grid g, double conductivity)
{
	vector <vector<double>> Hg(g.nodes.size());
	vector <vector<double>> Hl;

	for (int i = 0; i < Hg.size(); i++)
		Hg[i].resize(g.nodes.size());

	for (int i = 0; i < Hg.size(); i++)
		for (int j = 0; j < Hg[0].size(); j++)
			Hg[i][j] = 0.0;

	for (int k = 0; k < g.elements.size(); k++)
	{
		//ID globalne node'ow w elemencie
		//dla k = 0
		//	0	6	7	1
		//cout << g.elements[k].nodes[0]->ID << '\t';		//0
		//cout << g.elements[k].nodes[1]->ID << '\t';		//6
		//cout << g.elements[k].nodes[2]->ID << '\t';		//7
		//cout << g.elements[k].nodes[3]->ID << endl;		//1

		Hl = calculateMatrixHdV(g.elements[k], conductivity);

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Hg[g.elements[k].nodes[i]->ID][g.elements[k].nodes[j]->ID] += Hl[i][j];

			}
		}
	}
	return Hg;
}

vector<vector<double>> calculateGlobalMatrixC(Grid g, double c, double ro)
{
	vector <vector<double>> Cg(g.nodes.size());
	vector <vector<double>> Cl;

	for (int i = 0; i < Cg.size(); i++)
		Cg[i].resize(g.nodes.size());

	for (int i = 0; i < Cg.size(); i++)
		for (int j = 0; j < Cg[0].size(); j++)
			Cg[i][j] = 0.0;

	for (int k = 0; k < g.elements.size(); k++)
	{
		Cl = calculateMatrixC(g.elements[k], c, ro);

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Cg[g.elements[k].nodes[i]->ID][g.elements[k].nodes[j]->ID] += Cl[i][j];
			}
		}

	}

	return Cg;
}

vector<vector<double>> calculateGlobalMatrixHdS(Grid g, double convection, double heightGrid, double lenghGrid)
{
	vector <vector<double>> Hg(g.nodes.size());
	vector <vector<double>> Hl;

	for (int i = 0; i < Hg.size(); i++)
		Hg[i].resize(g.nodes.size());

	for (int i = 0; i < Hg.size(); i++)
		for (int j = 0; j < Hg[0].size(); j++)
			Hg[i][j] = 0.0;

	for (int k = 0; k < g.elements.size(); k++)
	{
		Hl = calculateMatrixHdS(g.elements[k], convection, heightGrid, lenghGrid);

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Hg[g.elements[k].nodes[i]->ID][g.elements[k].nodes[j]->ID] += Hl[i][j];
			}
		}
	}
	return Hg;
}

vector<double> calculateGlobalVectorP(Grid g, double convection, double heightGrid, double lenghGrid, double ambientTemperature)
{
	vector<double> Pg(g.nodes.size());
	vector<double> Pl;
	for (int i = 0; i < Pg.size(); i++)
		Pg[i] = 0.0;

	for (int k = 0; k < g.elements.size(); k++)
	{
		Pl = calculateVectorP(g.elements[k], convection, heightGrid, lenghGrid, ambientTemperature);
		for (int i = 0; i < 4; i++)
		{
			Pg[g.elements[k].nodes[i]->ID] += Pl[i];
		}
	}

	return Pg;
}