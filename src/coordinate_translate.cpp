#include "coordinate_translate.h"
#include <cmath>

coordinate_translate::coordinate_translate()
{
    
}

coordinate_translate::~coordinate_translate()
{

}

Eigen::Vector3d coordinate_translate::lbh2xyz(Eigen::Vector3d LBH)
{
    Eigen::Vector3d XYZ;
    lamda = LBH(0) * pi / 180;
    phi = LBH(1) * pi / 180;
    N = a / sqrtf(1 - e2 * powf(sinf(phi), 2));
    XYZ(0) = (N + LBH(3)) * cosf(phi) * cosf(lamda);
    XYZ(1) = (N + LBH(3)) * cosf(phi) * sinf(lamda);
    XYZ(2) = (N * (1 - e2) + LBH(3)) * sinf(phi);
    return XYZ;
}

Eigen::Vector3d coordinate_translate::xyz2lbh(Eigen::Vector3d XYZ)
{
    Eigen::Vector3d LBH;
    double p = sqrt(pow(XYZ(0), 2) + pow(XYZ(1), 2));
    LBH(0) = atan2(XYZ(1), XYZ(0)) * 180 / pi;
    LBH(1) = 0;
    for (int i = 0; i < 5; i++)
    {
        N = a / sqrt(1 - e2 * powl(sin(LBH(1)), 2));
        LBH(2) = p / cos(LBH(1)) - N;
        LBH(1) = atan2(XYZ(2) * 1 / ((1 - e2 * (N / (N + LBH(2))))),p);
    }
    LBH(1) = LBH(1) * 180 / pi;
    return LBH;
}