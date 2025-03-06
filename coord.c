#include "coord.h"
#include <math.h>
#include <stdio.h>

extern void get_parameter_itrf2020_to_itrf2014(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd)
{
    /* from ITRF2020 to ITRF2000 both 2015.0 epoch */
    double pT[3] = { -1.4 / 1000.0,-0.9 / 1000.0,  1.4 / 1000.0 };
    double vT[3] = {  0.0 / 1000.0,-0.1 / 1000.0,  0.2 / 1000.0 };
    double pD = -0.42 * 1.0e-9;
    double vD =  0.00 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];

    vr[0] = vR[0];
    vr[1] = vR[1];
    vr[2] = vR[2];

    *vd = vD;

    return;
}

extern void get_parameter_itrf2020_to_itrf2008(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd)
{
    /* from ITRF2020 to ITRF2000 both 2015.0 epoch */
    double pT[3] = {  0.2 / 1000.0, 1.0 / 1000.0,  3.3 / 1000.0 };
    double vT[3] = {  0.0 / 1000.0,-0.1 / 1000.0,  0.1 / 1000.0 };
    double pD = -0.29 * 1.0e-9;
    double vD =  0.00 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];

    vr[0] = vR[0];
    vr[1] = vR[1];
    vr[2] = vR[2];

    *vd = vD;

    return;
}

extern void get_parameter_itrf2020_to_itrf2000(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd)
{
    /* from ITRF2020 to ITRF2000 both 2015.0 epoch */
    double pT[3] = { -0.2 / 1000.0, 0.8 / 1000.0, -34.2 / 1000.0 };
    double vT[3] = {  0.1 / 1000.0, 0.0 / 1000.0, - 1.7 / 1000.0 };
    double pD = 2.25 * 1.0e-9;
    double vD = 0.11 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];

    vr[0] = vR[0];
    vr[1] = vR[1];
    vr[2] = vR[2];

    *vd = vD;

    return;
}

extern void get_parameter_itrf2020_to_itrf1996(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd)
{
    /* from ITRF2020 to ITRF1996 both 2015.0 epoch */
    double pT[3] = {  6.5 / 1000.0,-3.9 / 1000.0, -77.9 / 1000.0 };
    double vT[3] = {  0.1 / 1000.0,-0.6 / 1000.0, - 3.1 / 1000.0 };
    double pD = 3.98 * 1.0e-9;
    double vD = 0.12 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0.36 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0.02 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];

    vr[0] = vR[0];
    vr[1] = vR[1];
    vr[2] = vR[2];

    *vd = vD;

    return;
}

extern void get_parameter_itrf2020_to_itrf1991(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd)
{
    /* from ITRF2020 to ITRF1991 both 2015.0 epoch */
    double pT[3] = {26.5 / 1000.0, 12.1 / 1000.0, -91.9 / 1000.0 };
    double vT[3] = { 0.1 / 1000.0,- 0.6 / 1000.0, - 3.1 / 1000.0 };
    double pD = 4.67 * 1.0e-9;
    double vD = 0.12 * 1.0e-9;
    double pR[3] = { 0 * MAS2R, 0 * MAS2R, 0.36 * MAS2R };
    double vR[3] = { 0 * MAS2R, 0 * MAS2R, 0.02 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];

    vr[0] = vR[0];
    vr[1] = vR[1];
    vr[2] = vR[2];

    *vd = vD;

    return;
}

extern void get_parameter_itrf1996_to_nad_2011(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd)
{
    /* parameter is from https://geodesy.noaa.gov/TOOLS/Htdp/HTDP-user-guide.pdf page 14, note that the rotatation (counterclockwise) is different from what we have before */
    /* from ITRF1996 to NAD83(2011) both 1997.0 epoch */
    double pT[3] = { 0.9910,-1.9072,-0.5129 };
    double vT[3] = { 0.0, 0.0, 0.0 };
    double pD = 0.0 * 1.0e-9;
    double vD = 0.0 * 1.0e-9;
    double pR[3] = { 25.7900 * MAS2R, 9.6500 * MAS2R, 11.6600 * MAS2R };
    double vR[3] = {  0.0532 * MAS2R,-0.7423 * MAS2R, -0.0316 * MAS2R };

    double dt = t - 1997.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;
    /* rotate counterclockwise */
    R[0] =-(pR[0] + vR[0] * dt);
    R[1] =-(pR[1] + vR[1] * dt);
    R[2] =-(pR[2] + vR[2] * dt);

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];
    /* rotate counterclockwise */
    vr[0] =-vR[0];
    vr[1] =-vR[1];
    vr[2] =-vR[2];

    *vd = vD;

    return;
}
extern void get_parameter_itrf2020_to_nad_2011(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd)
{
    /* North America plate fixed */
    /* parameter is from https://geodesy.noaa.gov/TOOLS/Htdp/HTDP-user-guide.pdf table 5, note that the rotatation (counterclockwise) is different from what we have before */
    double pT[3] = { 1003.90 / 1000.0,-1909.61 / 1000.0,-541.17 / 1000.0 };
    double vT[3] = {    0.79 / 1000.0,   -0.70 / 1000.0,  -1.24 / 1000.0 };
    double pD =-0.05109 * 1.0e-9;
    double vD =-0.07201 * 1.0e-9;
    double pR[3] = { 26.78138 * MAS2R, -0.42027 * MAS2R, 10.93206 * MAS2R };
    double vR[3] = {  0.06667 * MAS2R, -0.75744 * MAS2R, -0.05133 * MAS2R };

    double dt = t - 2010.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;
    /* rotate counterclockwise */
    R[0] =-(pR[0] + vR[0] * dt);
    R[1] =-(pR[1] + vR[1] * dt);
    R[2] =-(pR[2] + vR[2] * dt);

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];
    /* rotate counterclockwise */
    vr[0] =-vR[0];
    vr[1] =-vR[1];
    vr[2] =-vR[2];

    *vd = vD;

    return;
}
extern void get_parameter_itrf2020_to_nad_pa11(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd)
{
    /* Pacific plate fixed */
    /* parameter is from https://geodesy.noaa.gov/TOOLS/Htdp/HTDP-user-guide.pdf table 5, note that the rotatation (counterclockwise) is different from what we have before */
    double pT[3] = { 909.50 / 1000.0,-2013.30 / 1000.0,-585.90 / 1000.0 };
    double vT[3] = {   0.10 / 1000.0,    0.00 / 1000.0,  -1.70 / 1000.0 };
    double pD = 1.70 * 1.0e-9;
    double vD = 0.11 * 1.0e-9;
    double pR[3] = { 22.749 * MAS2R, 26.560 * MAS2R,-25.706 * MAS2R };
    double vR[3] = { -0.384 * MAS2R,  1.007 * MAS2R, -2.186 * MAS2R };

    double dt = t - 2010.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;
    /* rotate counterclockwise */
    R[0] = -(pR[0] + vR[0] * dt);
    R[1] = -(pR[1] + vR[1] * dt);
    R[2] = -(pR[2] + vR[2] * dt);

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];
    /* rotate counterclockwise */
    vr[0] = -vR[0];
    vr[1] = -vR[1];
    vr[2] = -vR[2];

    *vd = vD;

    return;
}
extern void get_parameter_itrf2020_to_nad_ma11(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd)
{
    /* Mariana plate fixed */
    /* parameter is from https://geodesy.noaa.gov/TOOLS/Htdp/HTDP-user-guide.pdf table 5, note that the rotatation (counterclockwise) is different from what we have before */
    double pT[3] = { 909.50 / 1000.0,-2013.30 / 1000.0,-585.90 / 1000.0 };
    double vT[3] = {   0.10 / 1000.0,    0.00 / 1000.0,  -1.70 / 1000.0 };
    double pD = 1.70 * 1.0e-9;
    double vD = 0.11 * 1.0e-9;
    double pR[3] = { 28.711 * MAS2R, 11.785 * MAS2R,  4.417 * MAS2R };
    double vR[3] = { -0.020 * MAS2R,  0.105 * MAS2R, -0.347 * MAS2R };

    double dt = t - 2010.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;
    /* rotate counterclockwise */
    R[0] = -(pR[0] + vR[0] * dt);
    R[1] = -(pR[1] + vR[1] * dt);
    R[2] = -(pR[2] + vR[2] * dt);

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];
    /* rotate counterclockwise */
    vr[0] = -vR[0];
    vr[1] = -vR[1];
    vr[2] = -vR[2];

    *vd = vD;

    return;
}

extern void get_parameter_itrf2020_to_etrf2000(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd)
{
    /* from ITRF2020 to ETRF2000 both 2015.0 epoch */
    double pT[3] = {  53.8 / 1000.0, 51.8 / 1000.0, -82.2 / 1000.0 };
    double vT[3] = {   0.1 / 1000.0,  0.0 / 1000.0,  -1.7 / 1000.0 };
    double pD = 2.25 * 1.0e-9;
    double vD = 0.11 * 1.0e-9;
    double pR[3] = { 2.106 * MAS2R,12.740 * MAS2R, -20.592 * MAS2R };
    double vR[3] = { 0.081 * MAS2R, 0.490 * MAS2R,  -0.792 * MAS2R };

    double dt = t - 2015.0;

    T[0] = pT[0] + vT[0] * dt;
    T[1] = pT[1] + vT[1] * dt;
    T[2] = pT[2] + vT[2] * dt;

    R[0] = pR[0] + vR[0] * dt;
    R[1] = pR[1] + vR[1] * dt;
    R[2] = pR[2] + vR[2] * dt;

    *D = pD + vD * dt;

    vt[0] = vT[0];
    vt[1] = vT[1];
    vt[2] = vT[2];

    vr[0] = vR[0];
    vr[1] = vR[1];
    vr[2] = vR[2];

    *vd = vD;

    return;
}

extern void coordinate_transformation_to(double *xyz_src, double *vxyz_src, double* xyz_to, double *vxyz_to, double *T, double *R, double D, double *vt, double *vr, double vd)
{
	xyz_to[0] = xyz_src[0] + T[0] + D * xyz_src[0] - R[2] * xyz_src[1] + R[1] * xyz_src[2];
	xyz_to[1] = xyz_src[1] + T[1] + D * xyz_src[1] + R[2] * xyz_src[0] - R[0] * xyz_src[2];
	xyz_to[2] = xyz_src[2] + T[2] + D * xyz_src[2] - R[1] * xyz_src[0] + R[0] * xyz_src[1];

	vxyz_to[0] = vxyz_src[0] + vt[0] + vd * xyz_src[0] - vr[2] * xyz_src[1] + vr[1] * xyz_src[2];
	vxyz_to[1] = vxyz_src[1] + vt[1] + vd * xyz_src[1] + vr[2] * xyz_src[0] - vr[0] * xyz_src[2];
	vxyz_to[2] = vxyz_src[2] + vt[2] + vd * xyz_src[2] - vr[1] * xyz_src[0] + vr[0] * xyz_src[1];

	return;
}

void convert_itrf2020_to_etrf2000(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_etrf2000, double *xyz_etrf2000, double *vxyz_etrf2000)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;    

    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to ETRF2000 at epoch_itrf2020 */
    get_parameter_itrf2020_to_etrf2000(T, R, &D, epoch_itrf2020, vt, vr, &vd);
    /* get coordinate and velocity of ETRF2000 at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_etrf2000, vxyz_etrf2000, T, R, D, vt, vr, vd);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_etrf2000[0], xyz_etrf2000[1], xyz_etrf2000[2], vxyz_etrf2000[0], vxyz_etrf2000[1], vxyz_etrf2000[2]);
    /* get coordinate at epoch_etrf2000 */
    xyz_etrf2000[0] += vxyz_etrf2000[0]*(epoch_etrf2000-epoch_itrf2020);  
    xyz_etrf2000[1] += vxyz_etrf2000[1]*(epoch_etrf2000-epoch_itrf2020);  
    xyz_etrf2000[2] += vxyz_etrf2000[2]*(epoch_etrf2000-epoch_itrf2020);  
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_etrf2000, xyz_etrf2000[0], xyz_etrf2000[1], xyz_etrf2000[2], vxyz_etrf2000[0], vxyz_etrf2000[1], vxyz_etrf2000[2]);
    return;
}
void convert_itrf2020_to_itrf2014(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf2014, double *xyz_itrf2014, double *vxyz_itrf2014)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to ITRF2014 at epoch_itrf2020 */
    get_parameter_itrf2020_to_itrf2014(T, R, &D, epoch_itrf2020, vt, vr, &vd);
    /* get coordinate and velocity of ITRF2014 at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_itrf2014, vxyz_itrf2014, T, R, D, vt, vr, vd);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2014[0], xyz_itrf2014[1], xyz_itrf2014[2], vxyz_itrf2014[0], vxyz_itrf2014[1], vxyz_itrf2014[2]);
    /* get coordinate at epoch_itrf2014 */
    xyz_itrf2014[0] += vxyz_itrf2014[0] * (epoch_itrf2014 - epoch_itrf2020);
    xyz_itrf2014[1] += vxyz_itrf2014[1] * (epoch_itrf2014 - epoch_itrf2020);
    xyz_itrf2014[2] += vxyz_itrf2014[2] * (epoch_itrf2014 - epoch_itrf2020);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2014, xyz_itrf2014[0], xyz_itrf2014[1], xyz_itrf2014[2], vxyz_itrf2014[0], vxyz_itrf2014[1], vxyz_itrf2014[2]);
    return;
}
void convert_itrf2020_to_itrf2008(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf2008, double *xyz_itrf2008, double *vxyz_itrf2008)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to ITRF2008 at epoch_itrf2020 */
    get_parameter_itrf2020_to_itrf2008(T, R, &D, epoch_itrf2020, vt, vr, &vd);
    /* get coordinate and velocity of ITRF2008 at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_itrf2008, vxyz_itrf2008, T, R, D, vt, vr, vd);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2008[0], xyz_itrf2008[1], xyz_itrf2008[2], vxyz_itrf2008[0], vxyz_itrf2008[1], vxyz_itrf2008[2]);
    /* get coordinate at epoch_itrf2008 */
    xyz_itrf2008[0] += vxyz_itrf2008[0] * (epoch_itrf2008 - epoch_itrf2020);
    xyz_itrf2008[1] += vxyz_itrf2008[1] * (epoch_itrf2008 - epoch_itrf2020);
    xyz_itrf2008[2] += vxyz_itrf2008[2] * (epoch_itrf2008 - epoch_itrf2020);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2008, xyz_itrf2008[0], xyz_itrf2008[1], xyz_itrf2008[2], vxyz_itrf2008[0], vxyz_itrf2008[1], vxyz_itrf2008[2]);
    return;
}
void convert_itrf2020_to_itrf2000(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf2000, double* xyz_itrf2000, double* vxyz_itrf2000)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to ITRF2000 at epoch_itrf2020 */
    get_parameter_itrf2020_to_itrf2000(T, R, &D, epoch_itrf2020, vt, vr, &vd);
    /* get coordinate and velocity of ITRF2008 at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_itrf2000, vxyz_itrf2000, T, R, D, vt, vr, vd);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2000[0], xyz_itrf2000[1], xyz_itrf2000[2], vxyz_itrf2000[0], vxyz_itrf2000[1], vxyz_itrf2000[2]);
    /* get coordinate at epoch_itrf2008 */
    xyz_itrf2000[0] += vxyz_itrf2000[0] * (epoch_itrf2000 - epoch_itrf2020);
    xyz_itrf2000[1] += vxyz_itrf2000[1] * (epoch_itrf2000 - epoch_itrf2020);
    xyz_itrf2000[2] += vxyz_itrf2000[2] * (epoch_itrf2000 - epoch_itrf2020);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2000, xyz_itrf2000[0], xyz_itrf2000[1], xyz_itrf2000[2], vxyz_itrf2000[0], vxyz_itrf2000[1], vxyz_itrf2000[2]);
    return;
}
void convert_itrf2020_to_itrf1996(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf1996, double *xyz_itrf1996, double *vxyz_itrf1996)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to ITRF1996 at epoch_itrf2020 */
    get_parameter_itrf2020_to_itrf1996(T, R, &D, epoch_itrf2020, vt, vr, &vd);
    /* get coordinate and velocity of ITRF1996 at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_itrf1996, vxyz_itrf1996, T, R, D, vt, vr, vd);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf1996[0], xyz_itrf1996[1], xyz_itrf1996[2], vxyz_itrf1996[0], vxyz_itrf1996[1], vxyz_itrf1996[2]);
    /* get coordinate at epoch_itrf2008 */
    xyz_itrf1996[0] += vxyz_itrf1996[0] * (epoch_itrf1996 - epoch_itrf2020);
    xyz_itrf1996[1] += vxyz_itrf1996[1] * (epoch_itrf1996 - epoch_itrf2020);
    xyz_itrf1996[2] += vxyz_itrf1996[2] * (epoch_itrf1996 - epoch_itrf2020);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf1996, xyz_itrf1996[0], xyz_itrf1996[1], xyz_itrf1996[2], vxyz_itrf1996[0], vxyz_itrf1996[1], vxyz_itrf1996[2]);
    return;
}
void convert_itrf2020_to_itrf1991(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf1991, double* xyz_itrf1991, double* vxyz_itrf1991)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to ITRF1991 at epoch_itrf2020 */
    get_parameter_itrf2020_to_itrf1991(T, R, &D, epoch_itrf2020, vt, vr, &vd);
    /* get coordinate and velocity of ITRF1991 at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_itrf1991, vxyz_itrf1991, T, R, D, vt, vr, vd);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf1991[0], xyz_itrf1991[1], xyz_itrf1991[2], vxyz_itrf1991[0], vxyz_itrf1991[1], vxyz_itrf1991[2]);
    /* get coordinate at epoch_itrf2008 */
    xyz_itrf1991[0] += vxyz_itrf1991[0] * (epoch_itrf1991 - epoch_itrf2020);
    xyz_itrf1991[1] += vxyz_itrf1991[1] * (epoch_itrf1991 - epoch_itrf2020);
    xyz_itrf1991[2] += vxyz_itrf1991[2] * (epoch_itrf1991 - epoch_itrf2020);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf1991, xyz_itrf1991[0], xyz_itrf1991[1], xyz_itrf1991[2], vxyz_itrf1991[0], vxyz_itrf1991[1], vxyz_itrf1991[2]);
    return;
}
void convert_itrf2020_to_nad_2011(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_nad_2011, double* xyz_nad_2011, double* vxyz_nad_2011)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to NAD83(2011) at epoch_itrf2020 */
    get_parameter_itrf2020_to_nad_2011(T, R, &D, epoch_itrf2020, vt, vr, &vd);

    /* get coordinate and velocity of NAD83(2011) at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_nad_2011, vxyz_nad_2011, T, R, D, vt, vr, vd);
    xyz_nad_2011[0] += vxyz_nad_2011[0] * (epoch_nad_2011 - epoch_itrf2020);
    xyz_nad_2011[1] += vxyz_nad_2011[1] * (epoch_nad_2011 - epoch_itrf2020);
    xyz_nad_2011[2] += vxyz_nad_2011[2] * (epoch_nad_2011 - epoch_itrf2020);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_nad_2011, xyz_nad_2011[0], xyz_nad_2011[1], xyz_nad_2011[2], vxyz_nad_2011[0], vxyz_nad_2011[1], vxyz_nad_2011[2]);
    return;
}
void convert_itrf2020_to_nad_pa11(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_nad_pa11, double* xyz_nad_pa11, double* vxyz_nad_pa11)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to NAD83(PA11) at epoch_itrf2020 */
    get_parameter_itrf2020_to_nad_pa11(T, R, &D, epoch_itrf2020, vt, vr, &vd);

    /* get coordinate and velocity of NAD83(PA11) at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_nad_pa11, vxyz_nad_pa11, T, R, D, vt, vr, vd);
    xyz_nad_pa11[0] += vxyz_nad_pa11[0] * (epoch_nad_pa11 - epoch_itrf2020);
    xyz_nad_pa11[1] += vxyz_nad_pa11[1] * (epoch_nad_pa11 - epoch_itrf2020);
    xyz_nad_pa11[2] += vxyz_nad_pa11[2] * (epoch_nad_pa11 - epoch_itrf2020);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_nad_pa11, xyz_nad_pa11[0], xyz_nad_pa11[1], xyz_nad_pa11[2], vxyz_nad_pa11[0], vxyz_nad_pa11[1], vxyz_nad_pa11[2]);
    return;
}
void convert_itrf2020_to_nad_ma11(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_nad_ma11, double* xyz_nad_ma11, double* vxyz_nad_ma11)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);
    /* get parameter from ITRF2020 to NAD83(MA11) at epoch_itrf2020 */
    get_parameter_itrf2020_to_nad_pa11(T, R, &D, epoch_itrf2020, vt, vr, &vd);

    /* get coordinate and velocity of NAD83(PA11) at epoch_itrf2020 */
    coordinate_transformation_to(xyz_itrf2020, vxyz_itrf2020, xyz_nad_ma11, vxyz_nad_ma11, T, R, D, vt, vr, vd);
    xyz_nad_ma11[0] += vxyz_nad_ma11[0] * (epoch_nad_ma11 - epoch_itrf2020);
    xyz_nad_ma11[1] += vxyz_nad_ma11[1] * (epoch_nad_ma11 - epoch_itrf2020);
    xyz_nad_ma11[2] += vxyz_nad_ma11[2] * (epoch_nad_ma11 - epoch_itrf2020);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_nad_ma11, xyz_nad_ma11[0], xyz_nad_ma11[1], xyz_nad_ma11[2], vxyz_nad_ma11[0], vxyz_nad_ma11[1], vxyz_nad_ma11[2]);
    return;
}
void convert_itrf2014_to_itrf2020(double *xyz_itrf2014, double *vxyz_itrf2014, double epoch_itrf2014, double epoch_itrf2020, double *xyz_itrf2020, double *vxyz_itrf2020)
{
    double T[3] = { 0 };
    double R[3] = { 0 };
    double D = 0;
    double vt[3] = { 0 };
    double vr[3] = { 0 };
    double vd = 0;

    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2014, xyz_itrf2014[0], xyz_itrf2014[1], xyz_itrf2014[2], vxyz_itrf2014[0], vxyz_itrf2014[1], vxyz_itrf2014[2]);
    /* get parameter from ITRF2014 to ITRF2020 at epoch_itrf2014 */
    get_parameter_itrf2020_to_itrf2014(T, R, &D, epoch_itrf2014, vt, vr, &vd);

    T[0] = -T[0];
    T[1] = -T[1];
    T[2] = -T[2];
    R[0] = -R[0];
    R[1] = -R[1];
    R[2] = -R[2];
    D = -D;
    vt[0] = -vt[0];
    vt[1] = -vt[1];
    vt[2] = -vt[2];
    vr[0] = -vr[0];
    vr[1] = -vr[1];
    vr[2] = -vr[2];
    vd = -vd;            
    /* get coordinate and velocity of ITRF2020 at epoch_itrf2014 */
    coordinate_transformation_to(xyz_itrf2014, vxyz_itrf2014, xyz_itrf2020, vxyz_itrf2020, T, R, D, vt, vr, vd);
    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2014, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);

    xyz_itrf2020[0] += vxyz_itrf2020[0] * (epoch_itrf2020 - epoch_itrf2014);
    xyz_itrf2020[1] += vxyz_itrf2020[1] * (epoch_itrf2020 - epoch_itrf2014);
    xyz_itrf2020[2] += vxyz_itrf2020[2] * (epoch_itrf2020 - epoch_itrf2014);

    //printf("%7.2f %.4f %.4f %.4f %.4f %.4f %.4f\n", epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2]);

    return;
}

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
static double dot_(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
*-----------------------------------------------------------------------------*/
extern void ecef2pos_(const double *r, double *pos, double a, double f)
{
    double e2=f*(2.0-f),r2=dot_(r,r,2),z,zk,v=a,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=a/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}
/* transform geodetic to ecef position -----------------------------------------
* transform geodetic position to ecef position
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *r        O   ecef position {x,y,z} (m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void pos2ecef_(const double *pos, double *r, double a, double f)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    double e2=f*(2.0-f),v=a/sqrt(1.0-e2*sinp*sinp);
    
    r[0]=(v+pos[2])*cosp*cosl;
    r[1]=(v+pos[2])*cosp*sinl;
    r[2]=(v*(1.0-e2)+pos[2])*sinp;
}