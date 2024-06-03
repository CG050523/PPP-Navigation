#pragma once
#include <Eigen/Dense>

class coordinate_translate
{
public:
    double e2 = 0.00669437999013;
    double a = 6378137.0;
    double b = 6356752.3142;
    double pi = 3.1415926;
    double N = 0;
    double lamda = 0, phi = 0, h = 0;
    double x = 0, y = 0, z = 0;
    int press_flag = 0;
    coordinate_translate();
    ~coordinate_translate();
public:
    Eigen::Vector3d lbh2xyz(Eigen::Vector3d LBH);
    Eigen::Vector3d xyz2lbh(Eigen::Vector3d XYZ);
};
