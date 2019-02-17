#pragma once
#include <vector>
#include<iostream>
#include"Element.h"
#include"Node.h"

using namespace std;

class Grid
{
public:
	vector <Element> elements;
	vector <Node> nodes;

	Grid createGrid(int nH, int nL, double H, double L);
	void printGrid();

};