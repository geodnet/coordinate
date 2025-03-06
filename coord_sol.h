#pragma once
#include <string>
#include <map>
#include <algorithm>
#include <filesystem>
#include <cmath>
#include <cstring>

int is_na_station(std::string code);
int is_eu1_station(std::string code);
int is_eu2_station(std::string code);
void coordinate_difference(const double* xyz1, const double* xyz2, double* dxyz);
double baseline_length(const double* xyz1, const double* xyz2);
double lat2local(double lat, double* lat2north);
void BlhToCen(double lat, double lon, double c_en[3][3]);
void NED_to_XYZ(double lat, double lon, double* dned, double* dxyz);
void XYZ_to_NED(double lat, double lon, double* dxyz, double* dned);
void cov_convert(double lat, double lon, double* rms0, double* rms1);
int is_skip_year(int year);
int day_of_year(int year, int mon, int day);
double get_epoch(int year, int mon, int day);

int parse_fields(char* const buffer, char** val, char key);
std::string get_current_date_string();
std::string get_current_time_string();

int get_station_code(const char* name, double lat_deg, double lon_deg, std::string& code);
int get_station_vel(const char* name, double lat_deg, double lon_deg, double* vxyz);

int remove_lead(char* buffer);

#define MAXFIELD 100

struct coord_t
{
    std::string name;
    std::string code;
    std::string fname;
    std::string stime;
    std::string etime;
    std::string ctime;
    std::string type;
    std::string coord_name_itrf2020;
    std::string coord_name_regional;
    std::string coord_name_itrf2014;
    std::string coord_name_wgs84;
    std::string coord_name_itrf2020_2015;
    int flag;
    double rms;
    double amb_fix_rate;
    double epoch_itrf2020;
    double epoch_regional;
    double epoch_itrf2014;
    double epoch_wgs84;
    double epoch_itrf2020_2015;
    double xyz_itrf2020[3];
    double xyz_regional[3]; /* ecef xyz m */
    double xyz_itrf2014[3];
    double xyz_wgs84[3];
    double xyz_itrf2020_2015[3];
    double sigma95_xyz[3];
    double sigma95_NEU[3]; /* accuracy (95%) in NEU [m] */
    double blh_itrf2020[3];
    double blh_regional[3];
    double blh_itrf2014[3];
    double blh_wgs84[3];
    double blh_itrf2020_2015[3];
    double vxyz_itrf2020[3];
    double vxyz_regional[3];
    double vxyz_itrf2014[3];
    double vxyz_wgs84[3];
    coord_t();
    ~coord_t() {};
    coord_t(const coord_t& src);
    void reset();
    bool operator==(const coord_t& obj)
    {
        return name == obj.name && stime == obj.stime && type == obj.type;
    }
    bool operator< (const coord_t& obj)
    {
        if (name == obj.name)
            return name < obj.name;
        else if (stime == obj.stime)
        {
            return type < obj.type;
        }
        else
        {
            return stime < obj.stime;
        }
    }
    std::string get_date_time() const
    {
        std::string date_time = stime + '-' + name;
        if (date_time.length() > 24)
            date_time = date_time.substr(0, 24);
        return date_time;
    }
    int read_json_file(const char* fname, const char* stnname);
    int read_scsv_file(char* buffer);
    int convert_coord();
    int read_from_opus_file(const char* opusfname);
    int read_from_nrcan_file(const char* nrcanfname);
};