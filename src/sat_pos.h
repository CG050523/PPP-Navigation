#pragma once

//double A = 0;// 半长轴（3，4）
//double toe = 0; // 星历参考时刻（4，1）
//double t = 0; // 待计算指定时刻
//double delta_N = 0; // 平运动差（2，3）
//double M0 = 0; // 参考历元平近点角（2，4）
//double e = 0; // 偏心率（3，2）
//double omega_earth = 0; // 地球自转角速度
//double omega0 = .835906122525D + 00; // 本周初始历元的升交点赤经（4，3）
//double omega = .101794964119D + 01; // 近地点角距（5，3）
//double omega_DOT = -.837356307823D - 08; // 升交点赤经变化率（5，4）
//double i0 = .958243252552D + 00; // 轨道倾角（5，1）
//double IDOT = .503592405216D - 10; // 轨道倾角变化率，即delta_i0 / delta_t（6，1）
//double Cis = .614672899246D - 07; // 轨道倾角的正弦调和项改正的振幅（弧度）（4，4）
//double Cic = -.800937414169D - 07; // 轨道倾角的余弦调和项改正的振幅（弧度）（4，2）
//double Cus = .462494790554D - 05; // 纬度幅角的正弦调和项改正的振幅（弧度）（3，3）
//double Cuc = -.198744237423D - 05; // 纬度幅角的余弦调和项改正的振幅（弧度）（3，1）
//double Crs = -.405312500000D + 02; // 轨道半径的正弦调和项改正的振幅（m）（2，2）
//double Crc = .288406250000D + 03; // 轨道半径的余弦调和项改正的振幅（m）（6，2）
//double GM = 3.9860047e14; // 地球引力常数，设定值为 3.9860047D + 14
//
//% 计算平均运动角速度 n0
//n0 = sqrt(GM / (A * A * A));
//% 计算相对于星历参考历元的时间 tk（并做判断）
//tk = t - toe;
//if tk > 302400
//tk = tk - 608400;
//end
//if tk < -302400
//    tk = tk + 608400;
//end
//% 对平均运动角速度进行改正 n
//n = n0 + delta_N;
//% 计算平近点角 Mk
//Mk = M0 + n * tk;
//% 计算偏近点角 Ek
//syms Ek;
//[Ek] = vpasolve(Mk == Ek - e * sin(Ek));
//% 计算真近点角 f
//f = atan2((sqrt(1 - e * e) * sin(Ek)) / (1 - e * cos(Ek)), (cos(Ek) - e) / (1 - e * cos(Ek)));
//% 计算升交角距 u'
//u_discard = f + omega;
//% 计算升交角距改正数 σuk，向径改正数 σrk，轨道倾角改正数 σik
//sigma_uk = Cus * sin(2 * u_discard) + Cuc * cos(2 * u_discard);
//sigma_rk = Crs * sin(2 * u_discard) + Crc * cos(2 * u_discard);
//sigma_ik = Cis * sin(2 * u_discard) + Cic * cos(2 * u_discard);
//% 计算改正后的升交角距 uk，向径 rk 和轨道倾角 ik
//uk = u_discard + sigma_uk;
//rk = A * (1 - e * cos(Ek)) + sigma_rk;
//ik = i0 + sigma_ik + IDOT * tk;
//% 计算卫星在轨道平面上的位置(xk, yk)
//Xk = rk * cos(uk);
//Yk = rk * sin(uk);
//% 计算改正后的升交点经度 Lk
//Lk = omega0 + (omega_DOT - omega_earth) * t - omega_DOT * toe;
//% 计算卫星在地心坐标系下的坐标（X, Y, Z）
//X = Xk * cos(Lk) - Yk * cos(ik) * sin(Lk)
//Y = Xk * sin(Lk) + Yk * cos(ik) * cos(Lk)
//Z = Yk * sin(ik)