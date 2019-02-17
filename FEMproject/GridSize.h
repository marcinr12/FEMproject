#pragma once
#include<string>
#include<fstream>

class GridSize
{
public:
	double H = 0;
	double L = 0;
	int nH = 0;
	int nL = 0;

	bool getData();
};