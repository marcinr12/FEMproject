#include "jacobiMethod.h"

vector<vector<double>> podzielL(vector<vector<double>> tab)
{

	vector<vector<double>> l(tab.size());

	//uzupelnij zerami
	for (int i = 0; i < tab.size(); i++)
		for (int j = 0; j < tab.size(); j++)
			l[i].push_back(0);

	for (int i = 0; i < tab.size(); i++)
		for (int j = 0; j < i; j++)
			l[i][j] = tab[i][j];


	return l;
}

vector<vector<double>> podzielD(vector<vector<double>> tab)
{
	//stworz tablice
	vector<vector<double>> d(tab.size());

	//uzupelnij zerami
	for (int i = 0; i < d.size(); i++)
	{
		for (int j = 0; j < d.size(); j++)
			d[i].push_back(0);
	}

	for (int i = 0; i < tab.size(); i++)
		d[i][i] = tab[i][i];

	return d;
}

vector<vector<double>> podzielU(vector<vector<double>> tab)
{
	//stworz tablice
	vector<vector<double>> u(tab.size());

	//uzupelnij zerami
	for (int i = 0; i < tab.size(); i++)
		for (int j = 0; j < tab.size(); j++)
			u[i].push_back(0);

	for (int i = 0; i < tab.size(); i++)
		for (int j = 0; j < i; j++)
			u[j][i] = tab[j][i];


	return u;
}

vector<vector<double>> obliczN(vector<vector<double>> tab)
{
	//stworz tablice
	vector<vector<double>> n(tab.size());

	//uzupelnij zerami
	for (int i = 0; i < tab.size(); i++)
		for (int j = 0; j < tab.size(); j++)
			n[i].push_back(0);

	for (int i = 0; i < tab.size(); i++)
		n[i][i] = pow(tab[i][i], -1);

	return n;
}

vector<vector<double>> obliczSume(vector<vector<double>> tab1, vector<vector<double>> tab2)
{
	//stworz tablice
	vector<vector<double>> tab(tab1.size());

	//sumuj
	for (int i = 0; i < tab1.size(); i++)
		for (int j = 0; j < tab1.size(); j++)
			tab[i].push_back(tab1[i][j] + tab2[i][j]);

	return tab;
}

vector<vector<double>> obliczIloczyn(vector<vector<double>> tab1, vector<vector<double>> tab2)
{
	//stworz tablice
	vector<vector<double>> tab(tab1.size());

	double s;
	for (int i = 0; i < tab1.size(); i++)
		for (int j = 0; j < tab1.size(); j++)
		{
			s = 0;
			for (int k = 0; k < tab1.size(); k++)
				s += tab1[i][k] * tab2[k][j];
			tab[i].push_back(s);
		}

	return tab;
}

vector<vector<double>> obliczPrzeciwna(vector<vector<double>> tab)
{
	vector<vector<double>> u(tab.size());

	for (int i = 0; i < tab.size(); i++)
		for (int j = 0; j < tab.size(); j++)
			u[i].push_back(-1 * tab[i][j]);
	return u;
}

//		wsp => H	roz => {P}		rzad => int ile iteracji (dok³adnosc rozwiazan)
vector<double> metodaJacobiego(vector<vector<double>> wsp, vector<double> roz, vector <double> przyblizoneRozw, int rzad)
{
	//tworzenie i inicjalizacja wynikow
	vector<double> wyniki(wsp.size());
	for (int i = 0; i < wsp.size(); i++)
		wyniki[i] = przyblizoneRozw[i];

	//vector<double> tmp(wsp.size());
	//for (int i = 0; i < wsp.size(); i++)
	//	tmp[i] = przyblizoneRozw[i];

	vector<vector<double>> macierzD = podzielD(wsp);
	vector<vector<double>> macierzN = obliczN(macierzD);
	vector<vector<double>> macierzU = podzielU(wsp);
	vector<vector<double>>macierzL = podzielL(wsp);
	vector<vector<double>>macierz_N = obliczPrzeciwna(macierzN);

	vector<vector<double>> macierzSumaLU = obliczSume(macierzU, macierzL);
	vector<vector<double>>macierzM = obliczIloczyn(macierz_N, macierzSumaLU);


	for (int j = 0; j < rzad; j++)
	{
		for (int i = 0; i < wsp.size(); i++)
		{
			//tmp[i] = -1.0 * macierzN[i][i] * roz[i];
			wyniki[i] = -1.0 * macierzN[i][i] * roz[i];
			for (int j = 0; j < wsp.size(); j++)
			{
				//tmp[i] += macierzM[i][j] * wyniki[j];
				wyniki[i] += macierzM[i][j] * wyniki[j];
			}
		}

		/*for (int i = 0; i < wsp.size(); i++)
			wyniki[i] = tmp[i];*/
	}

	for (int j = 0; j < wsp.size(); j++)
	{
		wyniki[j] = -wyniki[j];
	}

	return wyniki;
}