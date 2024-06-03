#pragma once

#include <Eigen\Dense>
#include <fstream>
#include <vector>
#include <iostream>
#include "coordinate_translate.h"

using namespace Eigen;
using namespace std;

const double C_light = 299792458;
const double PI = 3.1415926;

class sat_data : coordinate_translate
{
private:
	typedef struct utc_time
	{
		int year = 0, mon = 0, day = 0, hour = 0, min = 0, sec = 0;
		int sat_num = 0;
	}utc_time;
	typedef struct delay_data
	{
		double sat_x = 0, sat_y = 0, sat_z = 0;
		double sat_ClockDelay = 0;//卫星钟差
		double relate_effect = 0;//相对论效应
		double earth_round = 0;//地球自转改正
		double tropospheric_delay = 0;//对流层改正
		double dis_fake = 0;//消电离层伪距
		double angle = 0;
	}delay_data;
	char sat_ID = 0;
	vector<utc_time> time_decode;
	vector<delay_data> delay_decode;

private:
	void data2Matrix(MatrixXd& matrix);
	auto readTxt(string filename, int end_line);
	double distance_calc(vector<delay_data>& data, int index, Vector3d& X);
	inline double deg2rad(double deg);
	void state_coeeficient(MatrixXd& A, Vector3d& X_0, VectorXd& L, int sat_num, int sat_num_all);
	Vector3d xyz2enu(Vector3d& xyz, Vector3d& station_XYZ);
public:
	void data_input(int rows, string file_location);
	int CountRows(string file_location);
	pair<MatrixXd, VectorXd> build_LS();
	pair<MatrixXd, VectorXd> build_RLS();
	pair<MatrixXd, VectorXd> kalman_filter();
};
