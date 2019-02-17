#include"GridSize.h"

using namespace std;

bool GridSize::getData()
{
	fstream plik("plik.txt", std::ios::in);
	if (plik.good())
	{
		string sH, sL, snH, snL;

		getline(plik, sH);
		getline(plik, sL);
		getline(plik, snH);
		getline(plik, snL);
		plik.close();

		this->H = stof(sH);
		this->L = stof(sL);
		this->nH = stoi(snH);
		this->nL = stoi(snL);
		return true;
	}
	return false;
}