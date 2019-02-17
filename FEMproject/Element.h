#pragma once
#include<vector>
#include "Node.h"

class Element
{
public:
	std::vector <Node*> nodes;	//wektor nodow

	Element();
	Element(Node *n1, Node *n2, Node *n3, Node *n4);
	//
	Element(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double c1, double c2, double c3, double c4);
	//
};