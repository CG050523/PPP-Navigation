#include <Eigen/Dense>
#include <iostream>
#include <ctime>
#include "Point_P.h"

using namespace Eigen;

string loc = "saterror_Ele.txt";
string out_loc = "out.txt";

int main()
{
	std::cerr << "reading... wait second\n";
	clock_t start_t, end_t;
	start_t = clock();
	sat_data sat;
	fstream coefficient;
	MatrixXd coefficient_matrix;
	sat.data_input(sat.CountRows(loc), loc);
	cerr.setf(ios::fixed);
	cerr.precision(5);
	std::cerr << "calculating...\n";
	auto result = sat.build_RLS();
	std::cerr << result.first << '\n' << result.second << '\n';
	auto kalman = sat.kalman_filter();
	std::cerr << kalman.first << '\n' << kalman.second << '\n';
	end_t = clock();
	std::cerr << "done!\ntakes: ";
	std::cerr << (end_t - start_t) * 1.0 / CLOCKS_PER_SEC << 's';
	return 0;
}