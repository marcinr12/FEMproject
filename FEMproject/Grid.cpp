#include"Grid.h"

Grid Grid::createGrid(int nH, int nL, double H, double L)
{
	int nodesAmount = nH * nL;
	int elementsAmount = (nH - 1) * (nL - 1);

	double x = 0, dx = 1.0 * H / (nH - 1);
	double y = 0, dy = 1.0 * L / (nL - 1);
	int nodeID = 0;
	for (int i = 0; i < nL; i++)	//uzupe³nianie siatki wêz³ami
	{
		for (int j = 0; j < nH; j++)
		{
			Node *n = new Node;
			n->ID = nodeID++;
			n->x = x;
			n->y = y;
			y += dy;
			this->nodes.push_back(*n);
		}
		y = 0;
		x += dx;
	}

	for (int i = 0; i < (nH - 1)*(nL - 1); ++i)
	{
		Element *e = new Element;
		this->elements.push_back(*e);
	}

	int n1 = 0;
	int m1 = nH;

	for (int k = 0; k < this->elements.size(); ++k)
	{
		if ((n1 + 1) % nH != 0)
		{
			this->elements[k].nodes.push_back(&(this->nodes[n1]));
			n1++;
			this->elements[k].nodes.push_back(&(this->nodes[m1]));
			m1++;
			this->elements[k].nodes.push_back(&(this->nodes[m1]));
			this->elements[k].nodes.push_back(&(this->nodes[n1]));
		}
		else
		{
			n1++;
			m1++;
			k--;
		}
	}

	return *this;
}

void Grid::printGrid()
{
	cout << "Printing grid..." << endl;
	for (int k = 0; k < this->elements.size(); ++k)
	{

		//cout << this->elements[k].nodes[0]->x << "\t" << this->elements[k].nodes[0]->y << endl;
		//cout << this->elements[k].nodes[1]->x << "\t" << this->elements[k].nodes[1]->y << endl;
		//cout << this->elements[k].nodes[2]->x << "\t" << this->elements[k].nodes[2]->y << endl;
		//cout << this->elements[k].nodes[3]->x << "\t" << this->elements[k].nodes[3]->y << endl;


		/*cout << this->elements[k].nodes[0]->ID << endl;
		cout << this->elements[k].nodes[1]->ID << endl;
		cout << this->elements[k].nodes[2]->ID << endl;
		cout << this->elements[k].nodes[3]->ID << endl;
		*/

		cout << this->elements[k].nodes[3]->ID << "\t";
		cout << this->elements[k].nodes[2]->ID << endl;
		cout << this->elements[k].nodes[0]->ID << "\t";
		cout << this->elements[k].nodes[1]->ID << endl;

		cout << endl << endl;
	}
}