#include "Element.h"

Element::Element() {}

Element::Element(Node *n1, Node *n2, Node *n3, Node *n4)
{
	nodes.push_back(n1);
	nodes.push_back(n2);
	nodes.push_back(n3);
	nodes.push_back(n4);
}

//
Element::Element(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double c1, double c2, double c3, double c4)
{
	Node *n1 = new Node(x1, y1, c1);
	Node *n2 = new Node(x2, y2, c2);
	Node *n3 = new Node(x3, y3, c3);
	Node *n4 = new Node(x4, y4, c4);

	nodes.push_back(n1);
	nodes.push_back(n2);
	nodes.push_back(n3);
	nodes.push_back(n4);
}
//