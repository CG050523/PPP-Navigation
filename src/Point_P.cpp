#include "Point_P.h"

auto sat_data::readTxt(string file_location, int end_line)
{
	ifstream text;
	text.open(file_location, ios::in);

	vector<string> strVec;
	for (int i = 0; i < end_line; i++)
	{
		string inbuf;
		getline(text, inbuf, '\n');
		strVec.push_back(inbuf);
	}
	text.close();
	return strVec;
}

int sat_data::CountRows(string file_location)
{
	int RowsOFcoefficient = 0;
	ifstream file;
	file.open(file_location, ios::in);
	string tmp;
	while (getline(file, tmp, '\n'))
	{
		RowsOFcoefficient++;
	}
	file.close();
	return RowsOFcoefficient;
}

void sat_data::data_input(int rows, string file_location)
{
	int start_line = 7;
	auto line_data = readTxt(file_location, rows);
	utc_time s_time;
	delay_data d_data;
	for (; start_line < rows; start_line++)
	{
		if (line_data[start_line].size() < 30)
		{
			s_time.year = atoi(data(line_data[start_line].substr(2, 4)));
			s_time.mon = atoi(data(line_data[start_line].substr(7, 2)));
			s_time.day = atoi(data(line_data[start_line].substr(10, 2)));
			s_time.hour = atoi(data(line_data[start_line].substr(14, 2)));
			s_time.min = atoi(data(line_data[start_line].substr(17, 2)));
			s_time.sec = atoi(data(line_data[start_line].substr(20, 2)));
			s_time.sat_num = atoi(data(line_data[start_line].substr(26, 2)));
			time_decode.push_back(s_time);
			continue;
		}
		d_data.sat_x = atof(data(line_data[start_line].substr(5, 14)));
		d_data.sat_y = atof(data(line_data[start_line].substr(23, 14)));
		d_data.sat_z = atof(data(line_data[start_line].substr(42, 14)));
		d_data.sat_ClockDelay = atof(data(line_data[start_line].substr(59, 14)));
		d_data.relate_effect = atof(data(line_data[start_line].substr(79, 15)));
		d_data.earth_round = atof(data(line_data[start_line].substr(102, 9)));
		d_data.tropospheric_delay = atof(data(line_data[start_line].substr(121, 7)));
		d_data.dis_fake = atof(data(line_data[start_line].substr(137, 14)));
		d_data.angle = atof(data(line_data[start_line].substr(166, 8)));
		delay_decode.push_back(d_data);
	}
}

void sat_data::data2Matrix(MatrixXd& matrix)
{
	for (int i = 0; i < delay_decode.size(); i++)
	{
		matrix(i, 0) = delay_decode[i].sat_x;
		matrix(i, 1) = delay_decode[i].sat_y;
		matrix(i, 2) = delay_decode[i].sat_z;
		matrix(i, 3) = delay_decode[i].sat_ClockDelay;
		matrix(i, 4) = delay_decode[i].relate_effect;
		matrix(i, 5) = delay_decode[i].earth_round;
		matrix(i, 6) = delay_decode[i].tropospheric_delay;
		matrix(i, 7) = delay_decode[i].dis_fake;
		matrix(i, 8) = delay_decode[i].angle;
	}
}

pair<MatrixXd, VectorXd> sat_data::build_LS()
{
	fstream xyz;
	xyz.setf(ios::fixed);
	xyz.precision(5);
	xyz.open("result_enu.txt", ios::out);
	pair<MatrixXd, VectorXd> out;
	Vector3d origin_xyz;
	origin_xyz << 4027881.341, 306998.792, 4919499.044;
	int sat_count = 0;//单个历元内观测到的卫星数
	int old_sat_count = 0;//已计算历元包含卫星数
	int next_T = static_cast<int>(time_decode.size());//欲计算历元数 1,10,20,30,40,50,100,200,300,400,500
	double sigma_0 = 0.3;
	MatrixXd P(10, 10);
	P = MatrixXd::Zero(10, 10);
	MatrixXd A;//系数值
	Matrix4d Qxx;
	Vector4d X_tmp = Vector4d::Zero();//用作判断的临时变量
	Vector4d X_T = Vector4d::Zero();//参数最终值
	Vector4d Delta_X_T = Vector4d::Zero();//参数改正数
	Vector3d X_0 = Vector3d::Zero();//参数估值
	VectorXd L;//常数项
	double dis = 0;//几何距离

	for (int index = 0; index < next_T; index++)
	{
		old_sat_count += time_decode[index].sat_num;
		sat_count = time_decode[index].sat_num;

		P.resize(sat_count, sat_count);
		P = MatrixXd::Zero(sat_count, sat_count);
		A.resize(sat_count, 4);
		A = MatrixXd::Zero(sat_count, 4);
		L.resize(sat_count);
		L = VectorXd::Zero(sat_count);
		Delta_X_T = Vector4d::Zero();
		X_tmp = Vector4d::Zero();
		X_T = Vector4d::Zero();
		//X_0 = Vector3d::Zero();
		//cerr << "%" << index << '\n';

		for (int i = 0; i < sat_count; i++)
		{
			P(i, i) = pow(sin(deg2rad(delay_decode[old_sat_count - sat_count + i].angle)), 2) / (sigma_0 * sigma_0);
		}
		while (true)
		{
			X_tmp = X_T;
			state_coeeficient(A, X_0, L, sat_count, old_sat_count);
			Qxx = (A.transpose() * P * A).inverse();
			Delta_X_T = Qxx * A.transpose() * P * L;
			X_T += Delta_X_T;
			if (sqrt(pow(Delta_X_T(0), 2) + pow(Delta_X_T(1), 2) + pow(Delta_X_T(2), 2)) < 1e-3) break;
			X_0 = X_T.head(3);
		}
		//auto enu = xyz2enu(X_T.head(3), origin_xyz);
		xyz << X_T.transpose() << '\n';
	}
	xyz.close();
	out.first = Qxx; out.second = X_T;
	return out;
}

pair<MatrixXd, VectorXd> sat_data::build_RLS()
{
	fstream xyz, xyz_enu;
	xyz.setf(ios::fixed);
	xyz.precision(5);
	xyz_enu.setf(ios::fixed);
	xyz_enu.precision(5);
	xyz.open("rls_output.txt", ios::out);
	xyz_enu.open("rls_enu_output.txt", ios::out);
	pair<MatrixXd, VectorXd> out;
	int sat_count = 0;//单个历元内观测到的卫星数
	int old_sat_count = 0;//已计算历元包含卫星数
	int next_T = static_cast<int>(time_decode.size());//欲计算历元数
	double sigma_0 = 0.3;
	MatrixXd P(10, 10);
	P = MatrixXd::Zero(10, 10);
	MatrixXd A(10, 4);//系数值
	Matrix3d Qxx;
	Matrix3d Qxx_tmp;
	Vector3d X_T = Vector3d::Zero();//参数最终值
	Vector3d X_0;
	X_0 << 4027895.82557, 307001.89377, 4919519.66484;//参数估值
	VectorXd L(10);//常数项

	Matrix3d B;
	Vector3d W;
	Matrix3d B_tmp;
	B_tmp = Matrix3d::Zero();
	Vector3d W_tmp;
	W_tmp = Vector3d::Zero();
	Matrix3d N_11;
	Vector3d N_12;
	Vector3d N_21;
	double N_22;
	Vector3d W_11;
	double W_22;

	double dis = 0;//几何距离
	old_sat_count = time_decode[0].sat_num;
	for (int index = 1; index < next_T; index++)
	{
		old_sat_count += time_decode[index].sat_num;
		sat_count = time_decode[index].sat_num;

		P.resize(sat_count, sat_count);
		P = MatrixXd::Zero(sat_count, sat_count);
		A.resize(sat_count, 4);
		A = MatrixXd::Zero(sat_count, 4);
		L.resize(sat_count);
		L = VectorXd::Zero(sat_count);
		X_T = X_0;
		for (int i = 0; i < sat_count; i++)
		{
			A(i, 3) = 1;
			P(i, i) = pow(sin(deg2rad(delay_decode[old_sat_count - sat_count + i].angle)), 2) / (sigma_0 * sigma_0);
		}
		state_coeeficient(A, X_0, L, sat_count, old_sat_count);
		N_11 = A.block(0, 0, sat_count, 3).transpose() * P * A.block(0, 0, sat_count, 3) + B_tmp;
		N_12 = A.block(0, 0, sat_count, 3).transpose() * P * A.block(0, 3, sat_count, 1);
		for (int j = 0; j < 3; j++)
		{
			N_21(j) = (A.block(0, 3, sat_count, 1).transpose() * P * A.block(0, 0, sat_count, 3))(0, j);
		}
		N_22 = (A.block(0, 3, sat_count, 1).transpose() * P * A.block(0, 3, sat_count, 1))(0, 0);
		W_11 = A.block(0, 0, sat_count, 3).transpose() * P * L + W_tmp;
		W_22 = (A.block(0, 3, sat_count, 1).transpose() * P * L)(0);
		B = N_11 - N_12 * (1 / N_22) * N_21.transpose();
		W = W_11 - N_12 * (1 / N_22) * W_22;
		X_T += B.inverse() * W;
		Qxx = B.inverse();
		B_tmp = B;
		W_tmp = W;
		auto enu = xyz2enu(X_T, X_0);
		xyz_enu << enu.transpose() << '\n';
		xyz << index + 1 << '\t' << X_T.transpose() << '\t' << sqrt(pow(X_T(0) - 4027881.341, 2) + pow(X_T(1) - 306998.792, 2) + pow(X_T(2) - 4919499.044, 2)) << '\n';
	}
	xyz.close();
	out.first = Qxx; out.second = X_T;
	return out;
}

Vector3d sat_data::xyz2enu(Vector3d& xyz, Vector3d& station_XYZ)
{
	Vector3d enu;
	auto station_LBH = xyz2lbh(station_XYZ);
	Matrix3d rotation;
	rotation.setIdentity();
	Vector3d xyz_origin;
	rotation(0, 0) = -sin(deg2rad(station_LBH(0)));
	rotation(0, 1) = cos(deg2rad(station_LBH(0)));
	rotation(1, 0) = -sin(deg2rad(station_LBH(1))) * cos(deg2rad(station_LBH(0)));
	rotation(1, 1) = -sin(deg2rad(station_LBH(1))) * sin(deg2rad(station_LBH(0)));
	rotation(1, 2) = cos(deg2rad(station_LBH(1)));
	rotation(2, 0) = cos(deg2rad(station_LBH(1))) * cos(deg2rad(station_LBH(0)));
	rotation(2, 1) = cos(deg2rad(station_LBH(1))) * sin(deg2rad(station_LBH(0)));
	rotation(2, 2) = sin(deg2rad(station_LBH(1)));
	xyz_origin(0) = xyz(0) - station_XYZ(0);
	xyz_origin(1) = xyz(1) - station_XYZ(1);
	xyz_origin(2) = xyz(2) - station_XYZ(2);
	//xyz_origin(3) = 1;
	enu = rotation * xyz_origin;
	return enu;
}

pair<MatrixXd, VectorXd> sat_data::kalman_filter()
{
	fstream xyz;
	xyz.setf(ios::fixed);
	xyz.precision(5);
	xyz.open("kalman_output.txt", ios::out);
	pair<MatrixXd, VectorXd> out;
	int sat_count = 0;//单个历元内观测到的卫星数
	int old_sat_count = 0;//已计算历元包含卫星数
	int next_T = static_cast<int>(time_decode.size());//欲计算历元数
	double sigma_0 = 0.3;
	MatrixXd P(10, 10);
	P = MatrixXd::Zero(10, 10);
	MatrixXd K(4, 10);
	MatrixXd A;//系数值
	Matrix4d Qxx;
	Matrix4d Qxx_tmp;
	Qxx_tmp.setIdentity();
	Qxx_tmp *= 30;
	Vector4d X_tmp;//用作判断的临时变量
	X_tmp << 4027881.341, 307001.89375, 4919499.044, 1398508.12732;
	Vector4d X_T = Vector4d::Zero();//参数最终值
	Vector4d Delta_X_T = Vector4d::Zero();//参数改正数
	Vector3d X_0;
	X_0 = Vector3d::Zero();
	VectorXd L;//常数项
	Matrix4d I;
	I.setIdentity();
	double dis = 0;//几何距离

	for (int index = 1; index < next_T; index++)
	{
		old_sat_count += time_decode[index].sat_num;
		sat_count = time_decode[index].sat_num;
		P.resize(sat_count, sat_count);
		P = MatrixXd::Zero(sat_count, sat_count);
		A.resize(sat_count, 4);
		A = MatrixXd::Zero(sat_count, 4);
		L.resize(sat_count);
		L = VectorXd::Zero(sat_count);
		//Delta_X_T = Vector4d::Zero();
		K.resize(4, sat_count);

		X_0 = X_tmp.head(3);
		for (int i = 0; i < sat_count; i++)
		{
			P(i, i) = pow(sin(deg2rad(delay_decode[old_sat_count - sat_count + i].angle)), 2) / (sigma_0 * sigma_0);
		}
		state_coeeficient(A, X_0, L, sat_count, old_sat_count);
		K = Qxx_tmp * A.transpose() * (A * Qxx_tmp * A.transpose() + P.inverse()).inverse();
		//cerr << K << '\n' << '\n' << L << '\n' << '\n';
		auto Qxx_test = (A.transpose() * P * A).inverse();
		Delta_X_T(3) = (Qxx_test * A.transpose() * P * L)(3);
		X_T = X_tmp + K * (L - A * Delta_X_T);
		Delta_X_T = X_T - X_tmp;
		//cerr << A * X_tmp << '\n' << '\n';
		Qxx = (I - K * A) * Qxx_tmp;
		//cerr << Qxx << '\n' << '\n';
		xyz << index + 1 << '\t' << X_T.head(3).transpose() << '\t' << sqrt(pow(X_T(0) - 4027881.341, 2) + pow(X_T(1) - 306998.792, 2) + pow(X_T(2) - 4919499.044, 2)) << '\n';
		X_tmp = X_T;
		Qxx_tmp = Qxx;
		Qxx_tmp(3, 3) = 30;
	}
	//cout << A << '\n';
	xyz.close();
	out.first = Qxx; out.second = X_T;
	return out;
}

void sat_data::state_coeeficient(MatrixXd& A, Vector3d& X_0, VectorXd& L, int sat_num, int sat_num_all)
{
	auto index = sat_num_all - sat_num;
	double dis = 0;
	for (int i = sat_num_all - sat_num; i < sat_num_all; i++)
	{
		dis = distance_calc(delay_decode, i, X_0);
		A(i - index, 0) = -(delay_decode[i].sat_x - X_0(0)) / dis;
		A(i - index, 1) = -(delay_decode[i].sat_y - X_0(1)) / dis;
		A(i - index, 2) = -(delay_decode[i].sat_z - X_0(2)) / dis;
		A(i - index, 3) = 1;
		L(i - index) = (delay_decode[i].dis_fake - dis -
			delay_decode[i].earth_round + delay_decode[i].sat_ClockDelay * C_light +
			delay_decode[i].relate_effect * C_light - delay_decode[i].tropospheric_delay);
	}
}

inline double sat_data::deg2rad(double deg)
{
	return deg * PI / 180;
}

double sat_data::distance_calc(vector<delay_data>& data, int index, Vector3d& X)
{
	double dis = sqrt(pow(data[index].sat_x - X(0), 2) + pow(data[index].sat_y - X(1), 2) +
		pow(data[index].sat_z - X(2), 2));
	return dis;
}



//暂时无用
auto CountCols(fstream& file, string location)
{
	int ColsOFcoefficient = 0;
	file.open(location, ios::in);
	string tmp;
	getline(file, tmp);
	vector <string> vec;
	int sublen = 0;
	int size = static_cast<int>(tmp.size());
	for (int i = 0; i < size; i++)
	{
		if (tmp[i] != ' ')
		{
			sublen++;
			continue;
		}
		vec.push_back(tmp.substr(i - sublen, sublen));
		sublen = 0;
	}
	if (sublen != 0)
	{
		vec.push_back(tmp.substr(size - sublen, sublen));
	}
	ColsOFcoefficient = static_cast<int> (vec.size());
	/*while (getline(file, tmp))
	{
		ColsOFcoefficient++;
	};*/
	file.close();
	return ColsOFcoefficient;
}