#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Grid.h"
using namespace std;

class Saver
{
public:
	static void save_to_file(string, vector<double>);
	static void save_to_file(string, int);
};