#include <iostream>
#include <random>
#include <fstream>
using namespace std;

constexpr int MIN = 0;
constexpr int MAX = 1;

constexpr int RAND_NUMS_TO_GENERATE = 2;

int main()
{
	random_device rd;
	default_random_engine eng(rd());
	uniform_int_distribution<int> distr(MIN, MAX);

	ofstream outpf;
	outpf.open ("outp.dat"); //opening output file

	if (!outpf) {
		cout << "No such output file" << endl;
	}
	else {
		cout << "Printing coordinates" << endl;
		for (int i = 0; i <= 39; ++i)
		{
			for (int j = 0; j <= 39; ++j) {
			int n = distr(eng);
			outpf << "	boxToCell { box (-6 $y_" << i << " $z_" << j << ") (64 $y_" << i+1 << " $z_" << j+1 << "); fieldValues (" << endl;
            if (n == 0) {
                outpf << "	volScalarFieldValue phi_alpha 1.0" << endl;
                outpf << "	volScalarFieldValue phi_beta 0.0 ); }" << endl; //printed in setFieldsDict format
            }
            if (n == 1) {
                outpf << "	volScalarFieldValue phi_alpha 0.0" << endl;
                outpf << "	volScalarFieldValue phi_beta 1.0 ); }" << endl; //printed in setFieldsDict format
            }
			}
		}
	}
	outpf.close();

	return 0;
}
