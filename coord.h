#ifndef _COORD_H_
#define _COORD_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PI
#define PI          3.1415926535897932  /* pi */
#endif

#ifndef D2R
#define D2R         (PI/180.0)          /* deg to rad */
#endif

#ifndef R2D
#define R2D         (180.0/PI)          /* rad to deg */
#endif

#ifndef RE_WGS84
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#endif

#ifndef FE_WGS84
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */
#endif 

#ifndef RE_GRS80
#define RE_GRS80    6378137.0           /* earth semimajor axis (WGS84) (m) */
#endif

#ifndef FE_GRS80
#define FE_GRS80    (1.0/298.257222101) /* earth flattening (WGS84) */
#endif

#ifndef MAS2R
#define MAS2R       (0.001/3600.0*D2R)
#endif

/* transformation parameters are from https://itrf.ign.fr/en/solutions/transformations */
void get_parameter_itrf2020_to_itrf2014(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd);
void get_parameter_itrf2020_to_itrf2008(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd);
void get_parameter_itrf2020_to_itrf2000(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd);
void get_parameter_itrf2020_to_itrf1997(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd);
void get_parameter_itrf2020_to_itrf1996(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd);
void get_parameter_itrf2020_to_itrf1994(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd);
void get_parameter_itrf2020_to_itrf1991(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd);
void get_parameter_itrf1996_to_nad_2011(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd);
void get_parameter_itrf2020_to_nad_2011(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd);
void get_parameter_itrf2020_to_nad_pa11(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd);
void get_parameter_itrf2020_to_nad_ma11(double* T, double* R, double* D, double t, double* vt, double* vr, double* vd);
void get_parameter_itrf2020_to_etrf2000(double* T, double* R, double *D, double t, double *vt, double *vr, double *vd);

/* transformation */
void coordinate_transformation_to(double* xyz_src, double* vxyz_src, double* xyz_to, double* vxyz_to, double* T, double* R, double D, double* vt, double* vr, double vd);

void convert_itrf2020_to_etrf2000(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_etrf2000, double *xyz_etrf2000, double *vxyz_etrf2000);
void convert_itrf2020_to_itrf2014(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf2014, double *xyz_itrf2014, double *vxyz_itrf2014);
void convert_itrf2020_to_itrf2008(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf2008, double *xyz_itrf2008, double *vxyz_itrf2008);
void convert_itrf2020_to_itrf2000(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf2000, double *xyz_itrf2000, double *vxyz_itrf2000);
void convert_itrf2020_to_itrf1997(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf1997, double *xyz_itrf1997, double *vxyz_itrf1997);
void convert_itrf2020_to_itrf1996(double *xyz_itrf2020, double *vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf1996, double *xyz_itrf1996, double *vxyz_itrf1996);
void convert_itrf2020_to_itrf1994(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf1994, double* xyz_itrf1994, double* vxyz_itrf1994);
void convert_itrf2020_to_itrf1991(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_itrf1991, double* xyz_itrf1991, double* vxyz_itrf1991);
void convert_itrf2020_to_nad_2011(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_nad_2011, double* xyz_nad_2011, double* vxyz_nad_2011);
void convert_itrf2020_to_nad_pa11(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_nad_pa11, double* xyz_nad_pa11, double* vxyz_nad_pa11);
void convert_itrf2020_to_nad_ma11(double* xyz_itrf2020, double* vxyz_itrf2020, double epoch_itrf2020, double epoch_nad_ma11, double* xyz_nad_ma11, double* vxyz_nad_ma11);
void convert_itrf2014_to_itrf2020(double* xyz_itrf2014, double* vxyz_itrf2014, double epoch_itrf2014, double epoch_itrf2020, double* xyz_itrf2020, double* vxyz_itrf2020);
void convert_nad_2011_to_itrf2020(double* xyz_nad_2011, double* vxyz_nad_2011, double epoch_nad_2011, double epoch_itrf2020, double* xyz_itrf2020, double* vxyz_itrf2020);
void convert_nad_pa11_to_itrf2020(double* xyz_nad_2011, double* vxyz_nad_2011, double epoch_nad_2011, double epoch_itrf2020, double* xyz_itrf2020, double* vxyz_itrf2020);
void convert_nad_ma11_to_itrf2020(double* xyz_nad_2011, double* vxyz_nad_2011, double epoch_nad_2011, double epoch_itrf2020, double* xyz_itrf2020, double* vxyz_itrf2020);


/* ecef xyz to latitude,longitude and height => ITRF need to use RE_GRS80 and FE_GRS80 */
void ecef2pos_(const double *r, double *pos, double a, double f);
void pos2ecef_(const double *pos, double *r, double a, double f);

void lat2ned(double lat, double* lat2north, double* lat2east);
void blh2cen(double lat, double lon, double* cen);
void ned2xyz(double lat, double lon, double* dned, double* dxyz);
void xyz2ned(double lat, double lon, double* dxyz, double* dned);

#ifdef __cplusplus
}
#endif


#endif
