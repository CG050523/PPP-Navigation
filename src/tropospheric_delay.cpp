#include "tropospheric_delay.h"

tropospheric_delay::tropospheric_delay()
{
}

tropospheric_delay::~tropospheric_delay()
{
}

double tropospheric_delay::calc_ComponentDelay(double dis_before_tropospheric_fix)
{
	ato_P = 1013.25 * powl((1 - 2.2557 * receiver_height * 1e-5), 5.2568);
	ato_T = 15.0 - 6.5 * receiver_height * 1e-3 + 273.15;
	ato_water_p = 6.108 * expf((17.15 * ato_T - 4684.0) / (ato_T - 38.45)) * h_rel;
	dry_componentDelay = (0.0022768 * ato_P /
		(1.0 - 0.00266 * cosf(2 * receiver_lamda_u) - 0.00028 * receiver_height * 1e-3)) *
		(1 / cosf(PI / 2 - deg2rad(E)));
	wet_componentDelay = 0.0022768 * (1255 / ato_T + 0.05) *
		ato_water_p * (1 / cosf(PI / 2 - deg2rad(E)));
	componentDelay = dry_componentDelay + wet_componentDelay;
	dis_after = dis_before + componentDelay;
	return dis_after;
}

inline double tropospheric_delay::deg2rad(double deg)
{
	return deg * PI / 180;
}