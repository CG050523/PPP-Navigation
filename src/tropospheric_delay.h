#pragma once
#include <cmath>

class tropospheric_delay
{
private:
	double PI = 3.1415926;
	double receiver_phi_u = 113.91168121;//接收机经度
	double receiver_lamda_u = 22.52138978;//接收机纬度
	double receiver_height = 23.36;//接收机高度

public:
	double E = 0;//俯仰角
	double A = 0;//方位角
	double componentDelay = 0;//对流层延迟
	double dry_componentDelay = 0;//对流层延迟干
	double wet_componentDelay = 0;//对流层延迟湿
	double dis_before = 0;//改正前伪距
	double dis_after = 0;//改正后伪距
	double ato_P = 0;//大气压
	double ato_T = 0;//大气绝对温度
	double ato_water_p = 0;//水汽压
	float h_rel = 0.7;//相对湿度

	tropospheric_delay();
	~tropospheric_delay();

	double deg2rad(double deg);

public:
	double calc_ComponentDelay(double dis_before_tropospheric_fix);
};