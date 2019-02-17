#include "Saver.h"
void Saver::save_to_file(string path, vector<double> t0)
{
	fstream file;
	file.open(path, ios::app);
	if (file.good() == false) return;
	file << endl << "T0" << endl;
	file << "{" << endl;
	for (int i = 0; i < t0.size(); ++i)
	{
		file << "t " << t0[i] << endl;
	}
	file << "}" << endl;
	file.close();
	return;
}
void Saver::save_to_file(string path, int iterations)
{
	fstream file;
	file.open(path, ios::app);
	if (file.good() == false) return;
	file << endl << "Iterations" << endl;
	file << "{" << endl;
	file << "i " << iterations << endl;
	file << "}" << endl;
	file.close();
	return;
}