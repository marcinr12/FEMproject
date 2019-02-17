#include"primaryCalculations.h"

vector <double> N1calc()
{
	double ksi1 = -1.0 / sqrt(3), eta1 = ksi1;
	double ksi2 = 1.0 / sqrt(3), eta2 = eta1;
	double ksi3 = ksi2, eta3 = ksi3;
	double ksi4 = -ksi2, eta4 = eta3;

	vector <double> N1;
	N1.push_back(0.25 * (1 - ksi1) * (1 - eta1));
	N1.push_back(0.25 * (1 - ksi2) * (1 - eta2));
	N1.push_back(0.25 * (1 - ksi3) * (1 - eta3));
	N1.push_back(0.25 * (1 - ksi4) * (1 - eta4));
	return N1;
}

vector <double> N2calc()
{
	double ksi1 = -1.0 / sqrt(3), eta1 = ksi1;
	double ksi2 = 1.0 / sqrt(3), eta2 = eta1;
	double ksi3 = ksi2, eta3 = ksi3;
	double ksi4 = -ksi2, eta4 = eta3;

	vector <double> N2;
	N2.push_back(0.25 * (1 + ksi1) * (1 - eta1));
	N2.push_back(0.25 * (1 + ksi2) * (1 - eta2));
	N2.push_back(0.25 * (1 + ksi3) * (1 - eta3));
	N2.push_back(0.25 * (1 + ksi4) * (1 - eta4));

	return N2;
}

vector <double> N3calc()
{
	double ksi1 = -1.0 / sqrt(3), eta1 = ksi1;
	double ksi2 = 1.0 / sqrt(3), eta2 = eta1;
	double ksi3 = ksi2, eta3 = ksi3;
	double ksi4 = -ksi2, eta4 = eta3;

	vector <double> N3;
	N3.push_back(0.25 * (1 + ksi1) * (1 + eta1));
	N3.push_back(0.25 * (1 + ksi2) * (1 + eta2));
	N3.push_back(0.25 * (1 + ksi3) * (1 + eta3));
	N3.push_back(0.25 * (1 + ksi4) * (1 + eta4));

	return N3;
}

vector <double> N4calc()
{
	double ksi1 = -1.0 / sqrt(3), eta1 = ksi1;
	double ksi2 = 1.0 / sqrt(3), eta2 = eta1;
	double ksi3 = ksi2, eta3 = ksi3;
	double ksi4 = -ksi2, eta4 = eta3;

	vector <double> N4;
	N4.push_back(0.25 * (1 - ksi1) * (1 + eta1));
	N4.push_back(0.25 * (1 - ksi2) * (1 + eta2));
	N4.push_back(0.25 * (1 - ksi3) * (1 + eta3));
	N4.push_back(0.25 * (1 - ksi4) * (1 + eta4));

	return N4;
}

vector <double> calculateNsForPC(double ksi, double eta)
{
	vector <double> n;

	n.push_back(0.25*(1 - ksi)*(1 - eta));
	n.push_back(0.25*(1 + ksi)*(1 - eta));
	n.push_back(0.25*(1 + ksi)*(1 + eta));
	n.push_back(0.25*(1 - ksi)*(1 + eta));

	return n;
}

vector <double> dN_dKsi1Calc()
{
	double eta1 = -1.0 / sqrt(3);
	double eta2 = eta1;
	double eta3 = 1.0 / sqrt(3);
	double eta4 = eta3;

	vector <double> dN_dKsi1;
	dN_dKsi1.push_back(-0.25 * (1 - eta1));
	dN_dKsi1.push_back(-0.25 * (1 - eta2));
	dN_dKsi1.push_back(-0.25 * (1 - eta3));
	dN_dKsi1.push_back(-0.25 * (1 - eta4));

	return dN_dKsi1;
}

vector <double> dN_dKsi2Calc()
{
	double eta1 = -1.0 / sqrt(3);
	double eta2 = eta1;
	double eta3 = 1.0 / sqrt(3);
	double eta4 = eta3;

	vector <double> dN_dKsi2;
	dN_dKsi2.push_back(0.25 * (1 - eta1));
	dN_dKsi2.push_back(0.25 * (1 - eta2));
	dN_dKsi2.push_back(0.25 * (1 - eta3));
	dN_dKsi2.push_back(0.25 * (1 - eta4));

	return dN_dKsi2;
}

vector <double> dN_dKsi3Calc()
{
	double eta1 = -1.0 / sqrt(3);
	double eta2 = eta1;
	double eta3 = 1.0 / sqrt(3);
	double eta4 = eta3;

	vector <double> dN_dKsi3;
	dN_dKsi3.push_back(0.25 * (1 + eta1));
	dN_dKsi3.push_back(0.25 * (1 + eta2));
	dN_dKsi3.push_back(0.25 * (1 + eta3));
	dN_dKsi3.push_back(0.25 * (1 + eta4));

	return dN_dKsi3;
}

vector <double> dN_dKsi4Calc()
{
	double eta1 = -1.0 / sqrt(3);
	double eta2 = eta1;
	double eta3 = 1.0 / sqrt(3);
	double eta4 = eta3;

	vector <double> dN_dKsi4;
	dN_dKsi4.push_back(-0.25 * (1 + eta1));
	dN_dKsi4.push_back(-0.25 * (1 + eta2));
	dN_dKsi4.push_back(-0.25 * (1 + eta3));
	dN_dKsi4.push_back(-0.25 * (1 + eta4));

	return dN_dKsi4;
}

vector <double> dN_dEta1Calc()
{
	double ksi1 = -1.0 / sqrt(3);
	double ksi2 = 1.0 / sqrt(3);
	double ksi3 = ksi2;
	double ksi4 = -ksi2;




	vector<double> dN_dEta1;
	dN_dEta1.push_back(-0.25 * (1 - ksi1));
	dN_dEta1.push_back(-0.25 * (1 - ksi2));
	dN_dEta1.push_back(-0.25 * (1 - ksi3));
	dN_dEta1.push_back(-0.25 * (1 - ksi4));

	return dN_dEta1;
}

vector <double> dN_dEta2Calc()
{
	double ksi1 = -1.0 / sqrt(3);
	double ksi2 = 1.0 / sqrt(3);
	double ksi3 = ksi2;
	double ksi4 = -ksi2;

	vector<double> dN_dEta2;

	dN_dEta2.push_back(-0.25 * (1 + ksi1));
	dN_dEta2.push_back(-0.25 * (1 + ksi2));
	dN_dEta2.push_back(-0.25 * (1 + ksi3));
	dN_dEta2.push_back(-0.25 * (1 + ksi4));

	return dN_dEta2;
}

vector <double> dN_dEta3Calc()
{
	double ksi1 = -1.0 / sqrt(3);
	double ksi2 = 1.0 / sqrt(3);
	double ksi3 = ksi2;
	double ksi4 = -ksi2;

	vector<double> dN_dEta3;

	dN_dEta3.push_back(0.25 * (1 + ksi1));
	dN_dEta3.push_back(0.25 * (1 + ksi2));
	dN_dEta3.push_back(0.25 * (1 + ksi3));
	dN_dEta3.push_back(0.25 * (1 + ksi4));

	return dN_dEta3;
}

vector <double> dN_dEta4Calc()
{
	double ksi1 = -1.0 / sqrt(3);
	double ksi2 = 1.0 / sqrt(3);
	double ksi3 = ksi2;
	double ksi4 = -ksi2;

	vector<double> dN_dEta4;
	dN_dEta4.push_back(0.25 * (1 - ksi1));
	dN_dEta4.push_back(0.25 * (1 - ksi2));
	dN_dEta4.push_back(0.25 * (1 - ksi3));
	dN_dEta4.push_back(0.25 * (1 - ksi4));

	return dN_dEta4;

}

vector <double> J_1_1Cal(Element e, vector<double> dN_dKsi1, vector<double> dN_dKsi2, vector<double> dN_dKsi3, vector<double> dN_dKsi4)
{

	vector <double> J_1_1;
	J_1_1.push_back(dN_dKsi1[0] * e.nodes[0]->x + dN_dKsi2[0] * e.nodes[1]->x + dN_dKsi3[0] * e.nodes[2]->x + dN_dKsi4[0] * e.nodes[3]->x);
	J_1_1.push_back(dN_dKsi1[1] * e.nodes[0]->x + dN_dKsi2[1] * e.nodes[1]->x + dN_dKsi3[1] * e.nodes[2]->x + dN_dKsi4[1] * e.nodes[3]->x);
	J_1_1.push_back(dN_dKsi1[2] * e.nodes[0]->x + dN_dKsi2[2] * e.nodes[1]->x + dN_dKsi3[2] * e.nodes[2]->x + dN_dKsi4[2] * e.nodes[3]->x);
	J_1_1.push_back(dN_dKsi1[3] * e.nodes[0]->x + dN_dKsi2[3] * e.nodes[1]->x + dN_dKsi3[3] * e.nodes[2]->x + dN_dKsi4[3] * e.nodes[3]->x);

	return J_1_1;
}

vector <double> J_1_2Cal(Element e, vector<double> dN_dKsi1, vector<double> dN_dKsi2, vector<double> dN_dKsi3, vector<double> dN_dKsi4)
{
	vector <double> J_1_2;

	J_1_2.push_back(dN_dKsi1[0] * e.nodes[0]->y + dN_dKsi2[0] * e.nodes[1]->y + dN_dKsi3[0] * e.nodes[2]->y + dN_dKsi4[0] * e.nodes[3]->y);
	J_1_2.push_back(dN_dKsi1[1] * e.nodes[0]->y + dN_dKsi2[1] * e.nodes[1]->y + dN_dKsi3[1] * e.nodes[2]->y + dN_dKsi4[1] * e.nodes[3]->y);
	J_1_2.push_back(dN_dKsi1[2] * e.nodes[0]->y + dN_dKsi2[2] * e.nodes[1]->y + dN_dKsi3[2] * e.nodes[2]->y + dN_dKsi4[2] * e.nodes[3]->y);
	J_1_2.push_back(dN_dKsi1[3] * e.nodes[0]->y + dN_dKsi2[3] * e.nodes[1]->y + dN_dKsi3[3] * e.nodes[2]->y + dN_dKsi4[3] * e.nodes[3]->y);

	return J_1_2;
}

vector <double> J_2_1Cal(Element e, vector<double> dN_dEta1, vector<double> dN_dEta2, vector<double> dN_dEta3, vector<double> dN_dEta4)
{
	vector <double> J_2_1;
	J_2_1.push_back(dN_dEta1[0] * e.nodes[0]->x + dN_dEta2[0] * e.nodes[1]->x + dN_dEta3[0] * e.nodes[2]->x + dN_dEta4[0] * e.nodes[3]->x);
	J_2_1.push_back(dN_dEta1[1] * e.nodes[0]->x + dN_dEta2[1] * e.nodes[1]->x + dN_dEta3[1] * e.nodes[2]->x + dN_dEta4[1] * e.nodes[3]->x);
	J_2_1.push_back(dN_dEta1[2] * e.nodes[0]->x + dN_dEta2[2] * e.nodes[1]->x + dN_dEta3[2] * e.nodes[2]->x + dN_dEta4[2] * e.nodes[3]->x);
	J_2_1.push_back(dN_dEta1[3] * e.nodes[0]->x + dN_dEta2[3] * e.nodes[1]->x + dN_dEta3[3] * e.nodes[2]->x + dN_dEta4[3] * e.nodes[3]->x);

	return J_2_1;
}

vector <double> J_2_2Cal(Element e, vector<double> dN_dEta1, vector<double> dN_dEta2, vector<double> dN_dEta3, vector<double> dN_dEta4)
{
	vector <double> J_2_2;
	J_2_2.push_back(dN_dEta1[0] * e.nodes[0]->y + dN_dEta2[0] * e.nodes[1]->y + dN_dEta3[0] * e.nodes[2]->y + dN_dEta4[0] * e.nodes[3]->y);
	J_2_2.push_back(dN_dEta1[1] * e.nodes[0]->y + dN_dEta2[1] * e.nodes[1]->y + dN_dEta3[1] * e.nodes[2]->y + dN_dEta4[1] * e.nodes[3]->y);
	J_2_2.push_back(dN_dEta1[2] * e.nodes[0]->y + dN_dEta2[2] * e.nodes[1]->y + dN_dEta3[2] * e.nodes[2]->y + dN_dEta4[2] * e.nodes[3]->y);
	J_2_2.push_back(dN_dEta1[3] * e.nodes[0]->y + dN_dEta2[3] * e.nodes[1]->y + dN_dEta3[3] * e.nodes[2]->y + dN_dEta4[3] * e.nodes[3]->y);

	return J_2_2;
}

vector <double> detJCal(vector <double> J_1_1, vector <double> J_1_2, vector <double> J_2_1, vector <double> J_2_2)
{
	vector <double> detJ;
	detJ.push_back(J_1_1[0] * J_2_2[0] - J_2_1[0] * J_1_2[0]);
	detJ.push_back(J_1_1[1] * J_2_2[1] - J_2_1[1] * J_1_2[1]);
	detJ.push_back(J_1_1[2] * J_2_2[2] - J_2_1[2] * J_1_2[2]);
	detJ.push_back(J_1_1[3] * J_2_2[3] - J_2_1[3] * J_1_2[3]);

	return detJ;
}

vector < double> J_1_1_1Cal(vector<double>J_2_2, vector <double> detJ)
{
	vector <double> J_1_1_1;

	J_1_1_1.push_back(J_2_2[0] / detJ[0]);
	J_1_1_1.push_back(J_2_2[1] / detJ[1]);
	J_1_1_1.push_back(J_2_2[2] / detJ[2]);
	J_1_1_1.push_back(J_2_2[3] / detJ[3]);

	return J_1_1_1;
}

vector < double> J_1_1_2Cal(vector<double>J_1_2, vector <double> detJ)
{
	vector <double> J_1_1_2;
	J_1_1_2.push_back(J_1_2[0] / detJ[0]); ///-J
	J_1_1_2.push_back(J_1_2[1] / detJ[1]);
	J_1_1_2.push_back(J_1_2[2] / detJ[2]);
	J_1_1_2.push_back(J_1_2[3] / detJ[3]);

	return J_1_1_2;
}

vector < double> J_1_2_1Cal(vector<double>J_2_1, vector <double> detJ)
{
	vector <double> J_1_2_1;
	J_1_2_1.push_back(J_2_1[0] / detJ[0]);
	J_1_2_1.push_back(J_2_1[1] / detJ[1]);
	J_1_2_1.push_back(J_2_1[2] / detJ[2]);
	J_1_2_1.push_back(J_2_1[3] / detJ[3]);

	return J_1_2_1;
}

vector < double> J_1_2_2Cal(vector<double>J_1_1, vector <double> detJ)
{
	vector <double> J_1_2_2;
	J_1_2_2.push_back(J_1_1[0] / detJ[0]);
	J_1_2_2.push_back(J_1_1[1] / detJ[1]);
	J_1_2_2.push_back(J_1_1[2] / detJ[2]);
	J_1_2_2.push_back(J_1_1[3] / detJ[3]);

	return J_1_2_2;
}

vector <double> pc_dN_d(double j1, double j2, double dN_dK1, double dN_dK2, double dN_dK3, double dN_dK4, double dN_dE1, double dN_dE2, double dN_dE3, double dN_dE4)
{
	/*skladnik1 = J_1_1_1[0] * dN_dKsi1[0] + J_1_1_2[0] * dN_dEta1[0];
	skladnik2 = J_1_1_1[0] * dN_dKsi2[0] + J_1_1_2[0] * dN_dEta2[0];
	skladnik3 = J_1_1_1[0] * dN_dKsi3[0] + J_1_1_2[0] * dN_dEta3[0];
	skladnik4 = J_1_1_1[0] * dN_dKsi4[0] + J_1_1_2[0] * dN_dEta4[0];
	double pc1xTab[4] = { skladnik1,skladnik2, skladnik3, skladnik4 };*/

	vector <double> dN_d;
	dN_d.push_back(j1 * dN_dK1 + j2 * dN_dE1);
	dN_d.push_back(j1 * dN_dK2 + j2 * dN_dE2);
	dN_d.push_back(j1 * dN_dK3 + j2 * dN_dE3);
	dN_d.push_back(j1 * dN_dK4 + j2 * dN_dE4);

	return dN_d;
}
