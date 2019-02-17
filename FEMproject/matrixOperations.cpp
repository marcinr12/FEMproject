#include"matrixOperations.h"

vector <vector <double>> multiplyVectorByTransparent(vector <double> v)
{
	vector<vector<double>> toReturn(v.size());
	for (int i = 0; i < v.size(); i++)
		toReturn[i].resize(v.size());

	for (int i = 0; i < toReturn.size(); i++)
	{
		for (int j = 0; j < toReturn[0].size(); j++)
		{
			toReturn[i][j] = v[i] * v[j];
		}
	}
	return toReturn;
}

vector <double> multiplyVectorByNumber(vector <double> v, double d)
{
	vector<double> toReturn(v.size());

	for (int i = 0; i < toReturn.size(); i++)
	{
		toReturn[i] = v[i] * d;

	}
	return toReturn;

}

vector <vector <double>> multiplyMatrixByNumber(vector <vector <double>> v, double d)
{
	vector<vector<double>> toReturn(v.size());
	for (int i = 0; i < v.size(); i++)
		toReturn[i].resize(v.size());

	for (int i = 0; i < toReturn.size(); i++)
	{
		for (int j = 0; j < toReturn[0].size(); j++)
		{
			toReturn[i][j] = v[i][j] * d;
		}
	}
	return toReturn;
}

vector <double> multiplyMatrixByVectorVector(vector<vector<double>> m, vector <double> u)
{
	vector<double> toReturn(m.size());
	for (int i = 0; i < toReturn.size(); i++)
		toReturn[i] = 0;

	for (int i = 0; i < m.size(); i++)
		for (int j = 0; j < m[0].size(); j++)
			toReturn[i] += m[i][j] * u[j];

	return toReturn;
}

vector <vector <double>> divideMatrix(vector <vector <double>> v, double d)
{
	vector<vector<double>> toReturn(v.size());
	for (int i = 0; i < v.size(); i++)
		toReturn[i].resize(v.size());

	for (int i = 0; i < toReturn.size(); i++)
	{
		for (int j = 0; j < toReturn[0].size(); j++)
		{
			toReturn[i][j] = v[i][j] / d;
		}
	}
	return toReturn;
}

vector <vector <double>> sumMatrix(vector <vector <double>> v, vector <vector <double>> u)
{
	vector<vector<double>> toReturn(v.size());
	for (int i = 0; i < v.size(); i++)
		toReturn[i].resize(v.size());

	for (int i = 0; i < toReturn.size(); i++)
	{
		for (int j = 0; j < toReturn[0].size(); j++)
		{
			toReturn[i][j] = v[i][j] + u[i][j];
		}
	}
	return toReturn;
}

vector <double> sumVector(vector <double> v, vector <double> u)
{
	vector <double> toReturn(v.size());

	for (int i = 0; i < v.size(); i++)
		toReturn[i] = v[i] + u[i];

	return toReturn;
}

void printMatrix(vector <vector <double>> v)
{
	for (int i = 0; i < v.size(); i++)
	{
		for (int j = 0; j < v[0].size(); j++)
		{
			cout << v[i][j] << '\t';
		}
		cout << endl;
	}
	cout << endl;
}

void printVector(vector <double> v)
{
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << '\t';
	cout << endl;
}