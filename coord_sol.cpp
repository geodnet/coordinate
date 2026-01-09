#include "coord_sol.h"
#include "coord.h"

int is_na_station(std::string code) /* north america */
{
    return (
        code.find("AIA") != std::string::npos ||
        code.find("ATG") != std::string::npos ||
        code.find("ABW") != std::string::npos ||
        code.find("BRB") != std::string::npos ||
        code.find("BLZ") != std::string::npos ||
        code.find("BMU") != std::string::npos ||
        code.find("BES") != std::string::npos ||
        code.find("VGB") != std::string::npos ||
        code.find("CAN") != std::string::npos ||
        code.find("CYM") != std::string::npos ||
        code.find("CRI") != std::string::npos ||
        code.find("CUB") != std::string::npos ||
        code.find("CUW") != std::string::npos ||
        code.find("DMA") != std::string::npos ||
        code.find("DOM") != std::string::npos ||
        code.find("SLV") != std::string::npos ||
        code.find("GRL") != std::string::npos ||
        code.find("GRD") != std::string::npos ||
        code.find("GLP") != std::string::npos ||
        code.find("GTM") != std::string::npos ||
        code.find("HTI") != std::string::npos ||
        code.find("HND") != std::string::npos ||
        code.find("JAM") != std::string::npos ||
        code.find("MTQ") != std::string::npos ||
        //code.find("MEX") != std::string::npos || => USE ITRF2014(2010.0)
        code.find("MSR") != std::string::npos ||
        code.find("ANT") != std::string::npos ||
        code.find("NIC") != std::string::npos ||
        code.find("PAN") != std::string::npos ||
        code.find("PRI") != std::string::npos ||
        code.find("BLM") != std::string::npos ||
        code.find("KNA") != std::string::npos ||
        code.find("LCA") != std::string::npos ||
        code.find("MAF") != std::string::npos ||
        code.find("SPM") != std::string::npos ||
        code.find("VCT") != std::string::npos ||
        code.find("SXM") != std::string::npos ||
        code.find("BHS") != std::string::npos ||
        code.find("TTO") != std::string::npos ||
        code.find("TCA") != std::string::npos ||
        code.find("USA") != std::string::npos ||
        code.find("VIR") != std::string::npos ||
        code.find("BLZ") != std::string::npos ||
        code.find("CRI") != std::string::npos ||
        code.find("COS") != std::string::npos
        );
}
int is_sa_station(std::string code) /* south america */
{
    return (
        code.find("ARG") != std::string::npos ||
        code.find("BOL") != std::string::npos ||
        code.find("BRA") != std::string::npos ||
        code.find("CHL") != std::string::npos ||
        code.find("COL") != std::string::npos ||
        code.find("ECU") != std::string::npos ||
        code.find("FLK") != std::string::npos ||
        code.find("GUF") != std::string::npos ||
        code.find("GUY") != std::string::npos ||
        code.find("PRY") != std::string::npos ||
        code.find("PER") != std::string::npos ||
        code.find("SUR") != std::string::npos ||
        code.find("URY") != std::string::npos ||
        code.find("VEN") != std::string::npos
        );
}


int is_eu1_station(std::string code)
{
    return (
        code.find("AUT") != std::string::npos ||
        code.find("BEL") != std::string::npos ||
        code.find("BGR") != std::string::npos ||
        code.find("HRV") != std::string::npos ||
        code.find("CYP") != std::string::npos ||
        code.find("CZE") != std::string::npos ||
        code.find("DNK") != std::string::npos ||
        code.find("EST") != std::string::npos ||
        code.find("FIN") != std::string::npos ||
        code.find("FRA") != std::string::npos ||
        code.find("DEU") != std::string::npos ||
        code.find("GRC") != std::string::npos ||
        code.find("HUN") != std::string::npos ||
        code.find("IRL") != std::string::npos ||
        code.find("ITA") != std::string::npos ||
        code.find("LVA") != std::string::npos ||
        code.find("LTU") != std::string::npos ||
        code.find("LUX") != std::string::npos ||
        code.find("MLT") != std::string::npos ||
        code.find("NLD") != std::string::npos ||
        code.find("POL") != std::string::npos ||
        code.find("PRT") != std::string::npos ||
        code.find("ROU") != std::string::npos ||
        code.find("SVK") != std::string::npos ||
        code.find("SVN") != std::string::npos ||
        code.find("ESP") != std::string::npos ||
        code.find("SWE") != std::string::npos
        );
}


int is_eu2_station(std::string code)
{
    return (
        code.find("ALB") != std::string::npos ||
        code.find("AND") != std::string::npos ||
        code.find("ARM") != std::string::npos ||
        code.find("BLR") != std::string::npos ||
        code.find("BIH") != std::string::npos ||
        code.find("FRO") != std::string::npos ||
        //code.find("GEO") != std::string::npos ||
        code.find("GIB") != std::string::npos ||
        code.find("ISL") != std::string::npos ||
        code.find("IMN") != std::string::npos ||
        code.find("XKX") != std::string::npos ||
        code.find("LIE") != std::string::npos ||
        code.find("MKD") != std::string::npos ||
        code.find("MDA") != std::string::npos ||
        code.find("MCO") != std::string::npos ||
        code.find("MNE") != std::string::npos ||
        code.find("NOR") != std::string::npos ||
        code.find("RUS") != std::string::npos ||
        code.find("SMR") != std::string::npos ||
        code.find("SRB") != std::string::npos ||
        code.find("CHE") != std::string::npos ||
        code.find("TUR") != std::string::npos ||
        code.find("UKR") != std::string::npos ||
        code.find("GBR") != std::string::npos ||
        code.find("VAT") != std::string::npos
        );
}
void coordinate_difference(const double* xyz1, const double* xyz2, double* dxyz)
{
    dxyz[0] = xyz2[0] - xyz1[0];
    dxyz[1] = xyz2[1] - xyz1[1];
    dxyz[2] = xyz2[2] - xyz1[2];
    return;
}
double baseline_length(const double* xyz1, const double* xyz2)
{
    double dxyz[3] = { 0 };
    coordinate_difference(xyz1, xyz2, dxyz);
    return sqrt(dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2]);
}

double lat2local(double lat, double* lat2north)
{
    double f_WGS84 = (1.0 / 298.257223563);
    double e2WGS84 = (2.0 * f_WGS84 - f_WGS84 * f_WGS84);
    double slat = sin(lat);
    double clat = cos(lat);
    double one_e2_slat2 = 1.0 - e2WGS84 * slat * slat;
    double Rn = 6378137.0 / sqrt(one_e2_slat2);
    double Rm = Rn * (1.0 - e2WGS84) / (one_e2_slat2);
    *lat2north = Rm;
    return Rn * clat;
}

void BlhToCen(double lat, double lon, double c_en[3][3])
{
    c_en[0][0] = -sin(lat) * cos(lon);
    c_en[1][0] = -sin(lat) * sin(lon);
    c_en[2][0] = cos(lat);
    c_en[0][1] = -sin(lon);
    c_en[1][1] = cos(lon);
    c_en[2][1] = 0.0;
    c_en[0][2] = -cos(lat) * cos(lon);
    c_en[1][2] = -cos(lat) * sin(lon);
    c_en[2][2] = -sin(lat);
}

void NED_to_XYZ(double lat, double lon, double* dned, double* dxyz)
{
    double C_en[3][3] = { 0 };
    BlhToCen(lat, lon, C_en);
    dxyz[0] = C_en[0][0] * dned[0] + C_en[0][1] * dned[1] + C_en[0][2] * dned[2];
    dxyz[1] = C_en[1][0] * dned[0] + C_en[1][1] * dned[1] + C_en[1][2] * dned[2];
    dxyz[2] = C_en[2][0] * dned[0] + C_en[2][1] * dned[1] + C_en[2][2] * dned[2];
    return;
}
void XYZ_to_NED(double lat, double lon, double* dxyz, double* dned)
{
    double C_en[3][3] = { 0 };
    BlhToCen(lat, lon, C_en);
    dned[0] = C_en[0][0] * dxyz[0] + C_en[1][0] * dxyz[1] + C_en[2][0] * dxyz[2];
    dned[1] = C_en[0][1] * dxyz[0] + C_en[1][1] * dxyz[1] + C_en[2][1] * dxyz[2];
    dned[2] = C_en[0][2] * dxyz[0] + C_en[1][2] * dxyz[1] + C_en[2][2] * dxyz[2];
    return;
}

void cov_convert(double lat, double lon, double* rms0, double* rms1)
{
    double C_en[3][3] = { 0 };
    BlhToCen(lat, lon, C_en);
    /* ned = C_en'*xyz => cov(ned) = C_en'*cov(xyz)*C_en */
    double Q[3][3] = { 0 };
    Q[0][0] = rms0[0] * rms0[0];
    Q[1][1] = rms0[1] * rms0[1];
    Q[2][2] = rms0[2] * rms0[2];
    double AQ[3][3] = { 0 };
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            double s = 0;
            for (int k = 0; k < 3; ++k)
            {
                s += C_en[k][i] * Q[k][j];
            }
            AQ[i][j] = s;
        }
    }
    double q[3][3] = { 0 };
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            double s = 0;
            for (int k = 0; k < 3; ++k)
            {
                s += AQ[i][k] * C_en[k][j];
            }
            q[i][j] = s;
        }
    }
    rms1[0] = sqrt(q[0][0]);
    rms1[1] = sqrt(q[1][1]);
    rms1[2] = sqrt(q[2][2]);
    return;
}


int is_skip_year(int year)
{
    return ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0);
}

int day_of_year(int year, int mon, int day) 
{
    int totalDay = day;
    int dayPerMon[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    for (int monIndex = 0; monIndex < mon - 1; ++monIndex)
        totalDay += dayPerMon[monIndex];
    if (mon > 2 && is_skip_year(year)) ++totalDay;
    return totalDay;
}

int doytotime(int year, int doy, int *mon, int *day) 
{
    int totalDay = doy;
    int dayPerMon[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    int maxday = 0;

    *mon = 0;
    *day = 0;

    for (int i=0;i<12;++i)
    {
        maxday = dayPerMon[i];
        if (i==1 && is_skip_year(year)) maxday++;
        if (doy<=maxday)
        {
            *mon = i+1;
            *day = doy;
            break;
        }
        doy -= dayPerMon[i];
        if (i==1 && is_skip_year(year)) doy--;
    }

    return (*mon>0) ? 1 : 0;
}

double get_epoch(int year, int mon, int day)
{
    int doy = day_of_year(year, mon, day);
    int tday = is_skip_year(year) ? 366 : 365;
    return year + ((double)(doy)) / tday;
}

int parse_fields(char* const buffer, char** val, char key)
{
    char* p, * q;
    int n = 0;

    /* parse fields */
    for (p = buffer; *p && n < MAXFIELD; p = q + 1) {
        if (p == NULL) break;
        if ((q = strchr(p, key)) || (q = strchr(p, '\n')) || (q = strchr(p, '\r'))) {
            val[n++] = p; *q = '\0';
        }
        else break;
    }
    if (p) val[n++] = p;
    return n;
}

std::string get_current_date_string()
{
    time_t now = time(0);
    tm* gmtm = gmtime(&now);
    char buffer[30] = { 0 };
    sprintf(buffer, "%04i-%02i-%02i", gmtm->tm_year + 1900, gmtm->tm_mon + 1, gmtm->tm_mday);
    return std::string(buffer);
}
std::string get_current_time_string()
{
    time_t now = time(0);
    tm* gmtm = gmtime(&now);
    char buffer[255] = { 0 };
    sprintf(buffer, "%04i-%02i-%02i-%02i-%02i-%02i", gmtm->tm_year + 1900, gmtm->tm_mon + 1, gmtm->tm_mday, gmtm->tm_hour, (int)gmtm->tm_min, (int)gmtm->tm_sec);
    return std::string(buffer);
}
/* get country code */
int get_station_code(const char* name, double lat_deg, double lon_deg, std::string& code)
{
    int ret = 0;
    char buffer[255] = { 0 };
    sprintf(buffer, "%s-code.txt", name);
    std::string code_fname = std::string(buffer);
    /* node get_code.js -station=G001 -lat=39 -lon=-121" */
    std::filesystem::remove(code_fname);

    if (lon_deg > 180.0) lon_deg -= 360.0;
    sprintf(buffer, "node get_code.js -station=%s -lat=%.9f -lon=%.9f > %s ", name, lat_deg, lon_deg, code_fname.c_str());
    //printf("%s\n", buffer);
    system(buffer);

    FILE* fOUT = fopen(code_fname.c_str(), "r");
    char* val[MAXFIELD];
    while (fOUT && !feof(fOUT) && fgets(buffer, sizeof(buffer), fOUT))
    {
        int num = parse_fields(buffer, val, ',');
        if (num < 2) continue;
        if (strstr(val[1], name))
        {
            if (strstr(val[0], "undefined"))
                ret = 0;
            else
                ret = 1;
            code = std::string(val[0]);
            break;
        }
    }
    if (fOUT) fclose(fOUT);
    if (!ret)
    {
        if (strstr(name, "NLD000"))
        {
            code = "NLD";
            ret = 2;
        }
        else {
            //printf("%s failed in code\n", name);
        }
    }
    std::filesystem::remove(code_fname);
    return ret;
}
/* get station velocity */
int get_station_vel(const char* name, double lat_deg, double lon_deg, double* vxyz)
{
    int ret = 0;
#ifndef _WIN32
    char buffer[255] = { 0 };
    sprintf(buffer, "%s-vel-inp.txt", name); std::string inp_fname = std::string(buffer);
    sprintf(buffer, "%s-vel-cfg.txt", name); std::string cfg_fname = std::string(buffer);
    sprintf(buffer, "%s-vel-out.txt", name); std::string out_fname = std::string(buffer);

    FILE* fINP = fopen(inp_fname.c_str(), "w");

    if (lon_deg > 180.0) lon_deg -= 360.0;

    if (fINP)
    {
        fprintf(fINP, "%14.9f,%14.9f,%s\n", lat_deg, -lon_deg, name);
        fclose(fINP);
    }

    FILE* fCFG = fopen(cfg_fname.c_str(), "w");
    if (fCFG)
    {
        fprintf(fCFG, "2\n%s\n24\n5\n%s\n0\n0\n", out_fname.c_str(), inp_fname.c_str());
        fclose(fCFG);
    }

    std::filesystem::remove(out_fname);

    sprintf(buffer, "htdp < %s >/dev/null", cfg_fname.c_str());
    system(buffer);

    FILE* fOUT = fopen(out_fname.c_str(), "r");

    while (fOUT && !feof(fOUT) && fgets(buffer, sizeof(buffer), fOUT))
    {
        if (strlen(buffer) < 85) continue;
        double vn = atof(buffer + 60);
        double ve = atof(buffer + 69);
        double vu = atof(buffer + 77);
        buffer[24] = '\0';
        remove_lead(buffer);
        if (!strstr(buffer, name)) continue;
        double vned[3] = { vn / 1000.0, ve / 1000.0, -vu / 1000.0 };
        NED_to_XYZ(lat_deg * D2R, lon_deg * D2R, vned, vxyz);
        ret++;
        break;
    }
    if (fOUT) fclose(fOUT);
    std::filesystem::remove(cfg_fname);
    std::filesystem::remove(out_fname);
    std::filesystem::remove(inp_fname);
#endif
    return ret;
}

int remove_lead(char* buffer)
{
    size_t nlen = 0;
    while ((nlen = strlen(buffer)) > 0)
    {
        if (buffer[0] == '\'') buffer[0] = ' ';
        if (buffer[0] == '\"') buffer[0] = ' ';
        if (buffer[0] == ' ')
        {
            std::rotate(buffer + 0, buffer + 1, buffer + nlen);
            if (nlen > 0) buffer[nlen - 1] = '\0';
        }
        else
        {
            break;
        }
    }
    while ((nlen = strlen(buffer)) > 0)
    {
        if (buffer[nlen - 1] == '\'' || buffer[nlen - 1] == '\"' || buffer[nlen - 1] == ',' || buffer[nlen - 1] == ' ')
        {
            buffer[nlen - 1] = '\0';
        }
        else
        {
            break;
        }
    }
    return 1;
}

coord_t::coord_t()
{
    flag = 0;
    rms = 0;
    amb_fix_rate = 0;
    epoch_itrf2020 = 0;
    epoch_regional = 0;
    epoch_itrf2014 = 0;
    epoch_wgs84 = 0;
    epoch_itrf2020_2015 = 0;
    xyz_itrf2020[0] = xyz_itrf2020[1] = xyz_itrf2020[2] = 0;
    xyz_regional[0] = xyz_regional[1] = xyz_regional[2] = 0;
    xyz_itrf2014[0] = xyz_itrf2014[1] = xyz_itrf2014[2] = 0;
    xyz_wgs84[0] = xyz_wgs84[1] = xyz_wgs84[2] = 0;
    xyz_itrf2020_2015[0] = xyz_itrf2020_2015[1] = xyz_itrf2020_2015[2] = 0;
    sigma95_xyz[0] = sigma95_xyz[1] = sigma95_xyz[2] = 0;
    sigma95_NEU[0] = sigma95_NEU[1] = sigma95_NEU[2] = 0;
    blh_itrf2020[0] = blh_itrf2020[1] = blh_itrf2020[2] = 0;
    blh_regional[0] = blh_regional[1] = blh_regional[2] = 0;
    blh_itrf2014[0] = blh_itrf2014[1] = blh_itrf2014[2] = 0;
    blh_wgs84[0] = blh_wgs84[1] = blh_wgs84[2] = 0;
    blh_itrf2020_2015[0] = blh_itrf2020_2015[1] = blh_itrf2020_2015[2] = 0;
    vxyz_itrf2020[0] = vxyz_itrf2020[1] = vxyz_itrf2020[2] = 0;
    vxyz_regional[0] = vxyz_regional[1] = vxyz_regional[2] = 0;
    vxyz_itrf2014[0] = vxyz_itrf2014[1] = vxyz_itrf2014[2] = 0;
    vxyz_wgs84[0] = vxyz_wgs84[1] = vxyz_wgs84[2] = 0;
    dN = dE = dU = 0;
}
coord_t::coord_t(const coord_t &src)
{
    name = src.name;
    code = src.code;
    fname = src.fname;
    stime = src.stime;
    etime = src.etime;
    ctime = src.ctime;
    type = src.type;
    coord_name_itrf2020 = src.coord_name_itrf2020;
    coord_name_regional = src.coord_name_regional;
    coord_name_itrf2014 = src.coord_name_itrf2014;
    coord_name_wgs84 = src.coord_name_wgs84;
    coord_name_itrf2020_2015 = src.coord_name_itrf2020_2015;
    flag = src.flag;
    rms = src.rms;
    amb_fix_rate = src.amb_fix_rate;
    epoch_itrf2020 = src.epoch_itrf2020;
    epoch_regional = src.epoch_regional;
    epoch_itrf2014 = src.epoch_itrf2014;
    epoch_wgs84 = src.epoch_wgs84;
    epoch_itrf2020_2015 = src.epoch_itrf2020_2015;
    memcpy(xyz_itrf2020, src.xyz_itrf2020, sizeof(xyz_itrf2020));
    memcpy(xyz_regional, src.xyz_regional, sizeof(xyz_regional));
    memcpy(xyz_itrf2014, src.xyz_itrf2014, sizeof(xyz_itrf2014));
    memcpy(xyz_wgs84, src.xyz_wgs84, sizeof(xyz_wgs84));
    memcpy(xyz_itrf2020_2015, src.xyz_itrf2020_2015, sizeof(xyz_itrf2020_2015));

    memcpy(sigma95_xyz, src.sigma95_xyz, sizeof(sigma95_xyz));
    memcpy(sigma95_NEU, src.sigma95_NEU, sizeof(sigma95_NEU));

    memcpy(blh_itrf2020, src.blh_itrf2020, sizeof(blh_itrf2020));
    memcpy(blh_regional, src.blh_regional, sizeof(blh_regional));
    memcpy(blh_itrf2014, src.blh_itrf2014, sizeof(blh_itrf2014));
    memcpy(blh_wgs84, src.blh_wgs84, sizeof(blh_wgs84));
    memcpy(blh_itrf2020_2015, src.blh_itrf2020_2015, sizeof(blh_itrf2020_2015));

    memcpy(vxyz_itrf2020, src.vxyz_itrf2020, sizeof(vxyz_itrf2020));
    memcpy(vxyz_regional, src.vxyz_regional, sizeof(vxyz_regional));
    memcpy(vxyz_itrf2014, src.vxyz_itrf2014, sizeof(vxyz_itrf2014));
    memcpy(vxyz_wgs84, src.vxyz_wgs84, sizeof(vxyz_wgs84));

    dN = src.dN;
    dE = src.dE;
    dU = src.dU;
}
void coord_t::reset()
{
    name.clear();
    code.clear();
    fname.clear();
    stime.clear();
    etime.clear();
    ctime.clear();
    type.clear();
    coord_name_itrf2020.clear();
    coord_name_regional.clear();
    coord_name_itrf2014.clear();
    coord_name_wgs84.clear();
    coord_name_itrf2020_2015.clear();
    flag = 0;
    rms = 0;
    amb_fix_rate = 0;
    epoch_itrf2020 = 0;
    epoch_regional = 0;
    epoch_itrf2014 = 0;
    epoch_wgs84 = 0;
    epoch_itrf2020_2015 = 0;
    xyz_itrf2020[0] = xyz_itrf2020[1] = xyz_itrf2020[2] = 0;
    xyz_regional[0] = xyz_regional[1] = xyz_regional[2] = 0;
    xyz_itrf2014[0] = xyz_itrf2014[1] = xyz_itrf2014[2] = 0;
    xyz_wgs84[0] = xyz_wgs84[1] = xyz_wgs84[2] = 0;
    xyz_itrf2020_2015[0] = xyz_itrf2020_2015[1] = xyz_itrf2020_2015[2] = 0;
    sigma95_xyz[0] = sigma95_xyz[1] = sigma95_xyz[2] = 0;
    sigma95_NEU[0] = sigma95_NEU[1] = sigma95_NEU[2] = 0;
    blh_itrf2020[0] = blh_itrf2020[1] = blh_itrf2020[2] = 0;
    blh_regional[0] = blh_regional[1] = blh_regional[2] = 0;
    blh_itrf2014[0] = blh_itrf2014[1] = blh_itrf2014[2] = 0;
    blh_wgs84[0] = blh_wgs84[1] = blh_wgs84[2] = 0;
    blh_itrf2020_2015[0] = blh_itrf2020_2015[1] = blh_itrf2020_2015[2] = 0;
    vxyz_itrf2020[0] = vxyz_itrf2020[1] = vxyz_itrf2020[2] = 0;
    vxyz_regional[0] = vxyz_regional[1] = vxyz_regional[2] = 0;
    vxyz_itrf2014[0] = vxyz_itrf2014[1] = vxyz_itrf2014[2] = 0;
    vxyz_wgs84[0] = vxyz_wgs84[1] = vxyz_wgs84[2] = 0;

    dN = dE = dU = 0;
}
int coord_t::read_json_file(const char* fname, const char* stnname)
{
    int ret = 0;
    FILE* fJSON = fopen(fname, "r");
    char buffer[512] = { 0 };
    int is_name_found = 0;
    int is_itrf2020_found = 0;
    int is_itrf2020_2015_found = 0;
    int is_wgs84_found = 0;
    int is_regional_found = 0;
    int is_x = 0;
    int is_y = 0;
    int is_z = 0;
    char* temp = NULL;
    char* temp1 = NULL;

    FILE* fOUT = fopen("coords.json", "a");

    while (fJSON && !feof(fJSON) && fgets(buffer, sizeof(buffer), fJSON))
    {
        if (fOUT) fprintf(fOUT, "%s", buffer);
        //printf("%s", buffer);
        if (temp = strchr(buffer, '\n')) temp[0] = '\0';
        remove_lead(buffer);
        if (is_name_found)
        {
            if (is_itrf2020_found == 1)
            {
                if ((temp = strstr(buffer, "x:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vx:"))
                        vxyz_itrf2020[0] = atof(temp1 + 1);
                    else
                        xyz_itrf2020[0] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "y:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vy:"))
                        vxyz_itrf2020[1] = atof(temp1 + 1);
                    else
                        xyz_itrf2020[1] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "z:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vz:"))
                        vxyz_itrf2020[2] = atof(temp1 + 1);
                    else
                        xyz_itrf2020[2] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "iar:")) && (temp1 = strstr(temp, ":")))
                {
                    amb_fix_rate = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "epoch:")) && (temp1 = strstr(temp, ":")))
                {
                    epoch_itrf2020 = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "sigma95X:")) && (temp1 = strstr(temp, ":")))
                {
                    sigma95_xyz[0] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "sigma95Y:")) && (temp1 = strstr(temp, ":")))
                {
                    sigma95_xyz[1] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "sigma95Z:")) && (temp1 = strstr(temp, ":")))
                {
                    sigma95_xyz[2] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "dateObserved:")) && (temp1 = strstr(temp, ":")))
                {
                    remove_lead(temp1 + 1);
                    char* temp11 = strchr(temp1 + 1, ' ');
                    if (temp11) temp11[0] = '\0';
                    stime = std::string(temp1 + 1);
                }
                if ((temp = strstr(buffer, "dateComputed:")) && (temp1 = strstr(temp, ":")))
                {
                    remove_lead(temp1 + 1);
                    char* temp11 = strchr(temp1 + 1, ' ');
                    if (temp11) temp11[0] = '\0';
                    ctime = std::string(temp1 + 1);
                }
                if ((temp = strstr(buffer, "name:")) && (temp1 = strstr(temp, ":")))
                {
                    remove_lead(temp1 + 1);
                    coord_name_itrf2020 = std::string(temp1 + 1);
                }
                if (strstr(buffer, "}"))
                {
                    is_itrf2020_found = 2;
                }
            }
            else if (!is_itrf2020_found && strstr(buffer, "itrf2020:"))
            {
                is_itrf2020_found = 1;
            }
            if (is_itrf2020_2015_found == 1)
            {
                if (!strstr(buffer, "vx:") && (temp = strstr(buffer, "x:")) && (temp1 = strstr(temp, ":")))
                {
                    xyz_itrf2020_2015[0] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "vx:")) && (temp1 = strstr(temp, ":")))
                {
                    vxyz_itrf2020[0] = atof(temp1 + 1);
                }
                if (!strstr(buffer, "vy:") && (temp = strstr(buffer, "y:")) && (temp1 = strstr(temp, ":")))
                {
                    xyz_itrf2020_2015[1] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "vy:")) && (temp1 = strstr(temp, ":")))
                {
                    vxyz_itrf2020[1] = atof(temp1 + 1);
                }
                if (!strstr(buffer, "vz:") && (temp = strstr(buffer, "z:")) && (temp1 = strstr(temp, ":")))
                {
                    xyz_itrf2020_2015[2] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "vz:")) && (temp1 = strstr(temp, ":")))
                {
                    vxyz_itrf2020[2] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "name:")) && (temp1 = strstr(temp, ":")))
                {
                    remove_lead(temp1 + 1);
                    coord_name_itrf2020_2015 = std::string(temp1 + 1);
                    std::size_t nloc = coord_name_itrf2020_2015.find_last_of('(');
                    if (nloc != std::string::npos)
                    {
                        epoch_itrf2020_2015 = atof(coord_name_itrf2020_2015.substr(nloc + 1).c_str());
                    }
                    else
                    {
                        epoch_itrf2020_2015 = atof(coord_name_itrf2020_2015.c_str());
                    }
                }
                if (strstr(buffer, "}"))
                {
                    is_itrf2020_2015_found = 2;
                }
            }
            else if (!is_itrf2020_2015_found && strstr(buffer, "itrf2015:"))
            {
                is_itrf2020_2015_found = 1;
            }
            if (is_wgs84_found == 1)
            {
                if ((temp = strstr(buffer, "x:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vx:"))
                        vxyz_wgs84[0] = atof(temp1 + 1);
                    else
                        xyz_wgs84[0] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "y:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vy:"))
                        vxyz_wgs84[1] = atof(temp1 + 1);
                    else
                        xyz_wgs84[1] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "z:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vz:"))
                        vxyz_wgs84[2] = atof(temp1 + 1);
                    else
                        xyz_wgs84[2] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "name:")) && (temp1 = strstr(temp, ":")))
                {
                    remove_lead(temp1 + 1);
                    coord_name_wgs84 = std::string(temp1 + 1);
                    std::size_t nloc = coord_name_wgs84.find_last_of('(');
                    if (nloc != std::string::npos)
                    {
                        epoch_wgs84 = atof(coord_name_wgs84.substr(nloc + 1).c_str());
                    }
                    else
                    {
                        epoch_wgs84 = atof(coord_name_wgs84.c_str());
                    }
                }
                if (strstr(buffer, "}"))
                {
                    is_wgs84_found = 2;
                }
            }
            else if (!is_wgs84_found && strstr(buffer, "wgs84:"))
            {
                is_wgs84_found = 1;
            }

            if (is_regional_found == 1)
            {
                if ((temp = strstr(buffer, "x:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vx:"))
                        vxyz_regional[0] = atof(temp1 + 1);
                    else
                        xyz_regional[0] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "y:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vy:"))
                        vxyz_regional[1] = atof(temp1 + 1);
                    else
                        xyz_regional[1] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "z:")) && (temp1 = strstr(temp, ":")))
                {
                    if (strstr(buffer, "vz:"))
                        vxyz_regional[2] = atof(temp1 + 1);
                    else
                        xyz_regional[2] = atof(temp1 + 1);
                }
                if ((temp = strstr(buffer, "name:")) && (temp1 = strstr(temp, ":")))
                {
                    remove_lead(temp1 + 1);
                    coord_name_regional = std::string(temp1 + 1);
                    std::size_t nloc = coord_name_regional.find_last_of('(');
                    if (nloc != std::string::npos)
                    {
                        epoch_regional = atof(coord_name_regional.substr(nloc + 1).c_str());
                    }
                    else
                    {
                        epoch_regional = atof(coord_name_regional.c_str());
                    }
                }
                if (strstr(buffer, "}"))
                {
                    is_regional_found = 2;
                }
            }
            else if (!is_regional_found && strstr(buffer, "regional:"))
            {
                is_regional_found = 1;
            }
            if (name.size() > 0 && is_itrf2020_found == 2 && is_itrf2020_2015_found == 2 && is_wgs84_found == 2 && is_regional_found == 2)
            {
                if (code.length()<3)
                {
                    ecef2pos_(xyz_itrf2020, blh_itrf2020, RE_GRS80, FE_GRS80);
                    get_station_code(name.c_str(), blh_itrf2020[0] * R2D, blh_itrf2020[1] * R2D, code);
                }
                double vel3D = sqrt(vxyz_itrf2020[0]*vxyz_itrf2020[0]+vxyz_itrf2020[1]*vxyz_itrf2020[1]+vxyz_itrf2020[2]*vxyz_itrf2020[2]);
                if (vel3D>0.5)
                {
                    vxyz_itrf2020[0] /= 1000.0;
                    vxyz_itrf2020[1] /= 1000.0;
                    vxyz_itrf2020[2] /= 1000.0;
                }
                double vel3D_r = sqrt(vxyz_regional[0]*vxyz_regional[0]+vxyz_regional[1]*vxyz_regional[1]+vxyz_regional[2]*vxyz_regional[2]);
                if (vel3D_r>0.5)
                {
                    vxyz_regional[0] /= 1000.0;
                    vxyz_regional[1] /= 1000.0;
                    vxyz_regional[2] /= 1000.0;
                }
                double vel3D_wgs84 = sqrt(vxyz_wgs84[0]*vxyz_wgs84[0]+vxyz_wgs84[1]*vxyz_wgs84[1]+vxyz_wgs84[2]*vxyz_wgs84[2]);
                if (vel3D_wgs84>0.5)
                {
                    vxyz_wgs84[0] /= 1000.0;
                    vxyz_wgs84[1] /= 1000.0;
                    vxyz_wgs84[2] /= 1000.0;
                }
                type = "NRCAN";
                if (name == std::string(stnname))
                {
                    ret = 1;
                    break;
                }
                reset();
                is_name_found = 0;
            }
        }
        else if ((temp = strstr(buffer, "name:")) && (temp1 = strstr(temp, ":")))
        {
            remove_lead(temp1 + 1);
            name = std::string(temp1 + 1);
            is_name_found = 1;
            is_itrf2020_found = 0;
            is_itrf2020_2015_found = 0;
            is_wgs84_found = 0;
            is_regional_found = 0;
        }
    }
    if (fJSON) fclose(fJSON);
    if (fOUT) fclose(fOUT);
    return ret;
}
int coord_t::read_scsv_file(char* buffer)
{
    reset();
    int ret = 0;
    char* val[MAXFIELD];
    char* temp = strchr(buffer, '\n');
    if (temp) temp[0] = '\0';
    temp = strchr(buffer, '\r');
    if (temp) temp[0] = '\0';
    int num = parse_fields(buffer, val, ',');
    int loc = 0;
    if (num > 14)
    {
        for (int i = 0; i < num; ++i)
            remove_lead(val[i]);
        loc = 0;
        type = std::string(val[loc++]);
        code = std::string(val[loc++]);
        name = std::string(val[loc++]);
        stime = std::string(val[loc++]);
        ctime = std::string(val[loc++]);
        coord_name_itrf2020 = std::string(val[loc++]);
        epoch_itrf2020 = atof(val[loc++]);
        xyz_itrf2020[0] = atof(val[loc++]);
        xyz_itrf2020[1] = atof(val[loc++]);
        xyz_itrf2020[2] = atof(val[loc++]);
        amb_fix_rate = atof(val[loc++]);
        sigma95_xyz[0] = atof(val[loc++]);
        sigma95_xyz[1] = atof(val[loc++]);
        sigma95_xyz[2] = atof(val[loc++]);
        // vxyz_itrf2020[0] = atof(val[loc++]);
        // vxyz_itrf2020[1] = atof(val[loc++]);
        // vxyz_itrf2020[2] = atof(val[loc++]);
        fname = std::string(val[loc++]);
        if (strstr(val[0], "OPUS") && num > 28)
        {
            coord_name_regional = std::string(val[loc++]);
            epoch_regional = atof(val[loc++]);
            xyz_regional[0] = atof(val[loc++]);
            xyz_regional[1] = atof(val[loc++]);
            xyz_regional[2] = atof(val[loc++]);
            coord_name_itrf2014 = std::string(val[loc++]);
            epoch_itrf2014 = atof(val[loc++]);
            xyz_itrf2014[0] = atof(val[loc++]);
            xyz_itrf2014[1] = atof(val[loc++]);
            xyz_itrf2014[2] = atof(val[loc++]);
            rms = atof(val[loc++]);
        }
        ret = 1;
    }
    return ret;
}
int coord_t::convert_coord()
{
    /* convert the other coordinate from ITRF2020 and VEL + country code */
    int ret = 0;

#ifdef _WIN32
    return ret;
#endif

    double rms3D = sqrt(sigma95_xyz[0] * sigma95_xyz[0] + sigma95_xyz[1] * sigma95_xyz[1] + sigma95_xyz[2] * sigma95_xyz[2]);
    if (rms3D < 1.0e-6 || name.length() == 0 || fabs(xyz_itrf2020[0]) < 0.00001 || fabs(xyz_itrf2020[1]) < 0.00001 || fabs(xyz_itrf2020[2]) < 0.00001 || fabs(epoch_itrf2020) < 0.00001 || coord_name_itrf2020.c_str() == 0)
    {
        return ret;
    }

    ecef2pos_(xyz_itrf2020, blh_itrf2020, RE_GRS80, FE_GRS80);

    if (code.length()<3)
    {
        get_station_code(name.c_str(), blh_itrf2020[0] * R2D, blh_itrf2020[1] * R2D, code);
    }

    if (code.length() < 3) return 0;

    int vel_flag = sqrt(vxyz_itrf2020[0] * vxyz_itrf2020[0] + vxyz_itrf2020[1] * vxyz_itrf2020[1] + vxyz_itrf2020[2] * vxyz_itrf2020[2]) > 0.0001;

    if (!vel_flag)
    {
        double vxyz[3] = { 0 };
        if (get_station_vel(name.c_str(), blh_itrf2020[0] * R2D, blh_itrf2020[1] * R2D, vxyz))
        {
            vxyz_itrf2020[0] = vxyz[0];
            vxyz_itrf2020[1] = vxyz[1];
            vxyz_itrf2020[2] = vxyz[2];
            vel_flag = sqrt(vxyz_itrf2020[0] * vxyz_itrf2020[0] + vxyz_itrf2020[1] * vxyz_itrf2020[1] + vxyz_itrf2020[2] * vxyz_itrf2020[2]) > 0.0001;
        }        
    }

    if (!vel_flag) return 0;

    char buffer[255] = { 0 };

    ecef2pos_(xyz_itrf2020, blh_itrf2020, RE_GRS80, FE_GRS80);

    double sigma95_NEU_old[3] = { sigma95_NEU[0], sigma95_NEU[1], sigma95_NEU[2] };
    double rms_xyz[3] = { sigma95_xyz[0] / 2.0, sigma95_xyz[1] / 2.0, sigma95_xyz[2] / 2.0 };
    double rms_ned[3] = { 0 };

    cov_convert(blh_itrf2020[0], blh_itrf2020[1], rms_xyz, rms_ned);
    sigma95_NEU[0] = rms_ned[0] * 2.0;
    sigma95_NEU[1] = rms_ned[1] * 2.0;
    sigma95_NEU[2] = rms_ned[2] * 2.0;

    /* WGS84(2139) */
    epoch_wgs84 = vel_flag ? floor(epoch_itrf2020) + 0.5 : epoch_itrf2020;
    sprintf(buffer, "WGS84(G2139)(%7.2f)", epoch_wgs84);
    coord_name_wgs84 = std::string(buffer);
    /* convert ITRF2020 to WGS84 */
    convert_itrf2020_to_itrf2014(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_wgs84, xyz_wgs84, vxyz_wgs84);

    /* ITRF2020(2015) */
    epoch_itrf2020_2015 = vel_flag ? 2015.0 : epoch_itrf2020;
    sprintf(buffer, "ITRF2020(%7.2f)", epoch_itrf2020_2015);
    coord_name_itrf2020_2015 = std::string(buffer);
    /* predict to ITRF2020 2015 epoch */
    xyz_itrf2020_2015[0] = xyz_itrf2020[0] + vxyz_itrf2020[0] * (epoch_itrf2020_2015 - epoch_itrf2020);
    xyz_itrf2020_2015[1] = xyz_itrf2020[1] + vxyz_itrf2020[1] * (epoch_itrf2020_2015 - epoch_itrf2020);
    xyz_itrf2020_2015[2] = xyz_itrf2020[2] + vxyz_itrf2020[2] * (epoch_itrf2020_2015 - epoch_itrf2020);
    /* backup the regional coordinate */
    double xyz_regional_old[3] = { xyz_regional[0], xyz_regional[1], xyz_regional[2] };
    std::string coord_name_regional_old = coord_name_regional;
    /* regional */
    if (code.find("IND") != std::string::npos)     /* ITRF2008(2005.0) */
    {
        /* output ITRF2008(2005.0) */
        epoch_regional = vel_flag ? 2005.0 : epoch_itrf2020;
        sprintf(buffer, "ITRF2008(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2008 */
        convert_itrf2020_to_itrf2008(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    /* Georgia Geodetic Datum
    The Georgia Geodetic Datum (GGD) is based and aligned with the International
    Terrestrial Reference System (ITRS). The ellipsoid associated with the datum is the
    Geodetic Reference System 1980 (GRS80).
    The relationship between the new datum and the ITRS is realized through the
    ITRF2008/IGS08 coordinates for 12 CORS stations at epoch 2011.353 computed by
    the Institute Geographique National (IGN) in France [IGN, 2011]. Hence, all points
    coordinated in terms
    */
    else if (code.find("GEO") != std::string::npos) /* Georgia Geodetic Datum = ITRF2008(2011.353) */
    {
        /* output ITRF2008(2011.353) */
        epoch_regional = vel_flag ? 2011.353 : epoch_itrf2020;
        sprintf(buffer, "ITRF2008(%8.3f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2008 */
        convert_itrf2020_to_itrf2008(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("EGY") != std::string::npos)     /* ITRF2008(2011.811) */
    {
        /* output ITRF2008(2011.811) */
        epoch_regional = vel_flag ? 2011.811 : epoch_itrf2020;
        sprintf(buffer, "ITRF2008(%8.3f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2008 */
        convert_itrf2020_to_itrf2008(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("NGA") != std::string::npos)     /* Nigerian Geocentric Datum (NGD2012)=ITRF2008(2012.0) */
    {
        /* output ITRF2008(2012.0) */
        epoch_regional = vel_flag ? 2012.0 : epoch_itrf2020;
        sprintf(buffer, "NGD2012(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2008 */
        convert_itrf2020_to_itrf2008(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("PHL") != std::string::npos) /* PGD2020=ITRF2014(2020.044) https://pagenet.namria.gov.ph/AGN/ServicesAndFees.aspx */
    {
        /* output ITRF2014(2020.044) */
        epoch_regional = vel_flag ? 2020.044 : epoch_itrf2020;
        sprintf(buffer, "PGD2020(%8.3f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2014 */
        convert_itrf2020_to_itrf2014(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("MEX") != std::string::npos) /* ITRF2014(2010) */
    {
        /* output ITRF2014(2010.0) */
        epoch_regional = vel_flag ? 2010.0 : epoch_itrf2020;
        sprintf(buffer, "ITRF2014(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2014 */
        convert_itrf2020_to_itrf2014(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("KEN") != std::string::npos) /* ITRF2014 current epoch */
    {
        /* output ITRF2014 */
        epoch_regional = vel_flag ? epoch_itrf2020 : epoch_itrf2020;
        sprintf(buffer, "ITRF2014(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2014 */
        convert_itrf2020_to_itrf2014(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("AUS") != std::string::npos) /* GDA2020=ITRF2014(2020.0) */
    {
        /* output ITRF2014(2020.0) */
        epoch_regional = vel_flag ? 2020.0 : epoch_itrf2020;
        sprintf(buffer, "GDA2020(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2014 */
        convert_itrf2020_to_itrf2014(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("NZL") != std::string::npos) /* NZGD2000=ITRF96(2000.0) */
    {
        /* output ITRF96(2000.0) */
        epoch_regional = vel_flag ? 2000.0 : epoch_itrf2020;
        sprintf(buffer, "NZGD2000(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF96 */
        convert_itrf2020_to_itrf1996(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("TUR") != std::string::npos) /* TUREF ITRF96(2005.0) */
    {
        /* output ITRF96(2005.0) */
        epoch_regional = vel_flag ? 2005.0 : epoch_itrf2020;
        sprintf(buffer, "ITRF96(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF96 */
        convert_itrf2020_to_itrf1996(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("CHN") != std::string::npos) /* CGCS2000 ITRF97(2000.0) */
    {
        /* output ITRF97(2000.0) */
        epoch_regional = vel_flag ? 2000.0 : epoch_itrf2020;
        sprintf(buffer, "CGCS2000(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF97 */
        convert_itrf2020_to_itrf1996(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
#if 0
    else if (code.find("JPN") != std::string::npos) /* JGD2000=ITRF94(1997.0) */
    {
        /* output ITRF94(1997.0) */
        epoch_regional = vel_flag ? 1997.0 : epoch_itrf2020;
        sprintf(buffer, "JGD2000(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF94 */
        convert_itrf2020_to_itrf1994(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
#endif
    else if (code.find("JPN") != std::string::npos) /* JGD2011=ITRF2008 on 05/24/2011=DOY(144) */
    {
        /* output ITRF2008 */
        int doy = day_of_year(2011, 05, 24);
        epoch_regional = vel_flag ? 2011.0 + doy / 365.0 : epoch_itrf2020;
        sprintf(buffer, "JGD2011(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2008 */
        convert_itrf2020_to_itrf2008(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("IDN") != std::string::npos) /* IGRS2013=ITRF2008(2012.0) */
    { /* https://un-ggim-ap.org/sites/default/files/media/meetings/Plenary08/WG1_S2B_2%20Arief%20Syafii_Indonesia_Geodetic%20Reference%20Frame.pdf */
        /* output ITRF2008(2012.0) */
        int doy = day_of_year(2011, 05, 24);
        epoch_regional = vel_flag ? 2012.0 : epoch_itrf2020;
        sprintf(buffer, "IGRS2013(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2008 */
        convert_itrf2020_to_itrf2008(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }    
    else if (code.find("GUM") != std::string::npos) /* NAD83(MA11)(2010.0) */
    {
        /* output NAD83(MA11)(2010.0) */
        epoch_regional = vel_flag ? 2010.0 : epoch_itrf2020;
        sprintf(buffer, "NAD83(MA11)(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to NAD83(MA11)(2010.0) */
        convert_itrf2020_to_nad_ma11(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("CAN") != std::string::npos) /* NAD83(CSRS)(2010.0)  */
    {
        /* output NAD83(CSRS)(2010.0) */
        epoch_regional = vel_flag ? 2010.0 : epoch_itrf2020;
        sprintf(buffer, "NAD83(CSRS)(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to NAD83(CSRS)(2010.0) */
        convert_itrf2020_to_nad_2011(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("USA") != std::string::npos) /* NAD83(PA11)(2010.0) or NAD83(2011)(2010.0)*/
    {
        /* check HAWII or not */
        double xyz_hi[3] = { -5455305.454,	-2402063.624,	2262323.603 };
        double dxyz1[3] = { xyz_hi[0] - xyz_itrf2020[0], xyz_hi[1] - xyz_itrf2020[1], xyz_hi[2] - xyz_itrf2020[2] };
        double dist1 = sqrt(dxyz1[0] * dxyz1[0] + dxyz1[1] * dxyz1[1] + dxyz1[2] * dxyz1[2]) / 1000.0;
        if (dist1 < 650.0) /* within G31 650km range => HAWII */
        {
            /* output NAD83(PA11)(2010.0) */
            epoch_regional = vel_flag ? 2010.0 : epoch_itrf2020;
            sprintf(buffer, "NAD83(PA11)(%7.2f)", epoch_regional);
            coord_name_regional = std::string(buffer);
            /* convert ITRF2020 to NAD83(PA11)(2010.0) */
            convert_itrf2020_to_nad_pa11(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
        }
        else
        {
            /* output NAD83(2011)(2010.0) */
            epoch_regional = vel_flag ? 2010.0 : epoch_itrf2020;
            sprintf(buffer, "NAD83(2011)(%7.2f)", epoch_regional);
            coord_name_regional = std::string(buffer);
            /* convert ITRF2020 to NAD83(2011)(2010.0) */
            convert_itrf2020_to_nad_2011(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
        }
    }
    else if (code.find("ZAF") != std::string::npos)     /* ITRF1991(1994.0) */
    {
        /* output ITRF1991(1994.0) */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 1994;
        sprintf(buffer, "ITRF1991(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF1991 */
        convert_itrf2020_to_itrf1991(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("LKA") != std::string::npos)     /* WGS84(G730)=ITRF1991(1994.0) */
    {
        /* it seems that the WGS84 do not match the coordinate good, ask the LKA CORS for 4 stations' rinex data in 2024.026
        *  and get the ITRF2020 coordinate, then compare their official coordinate, and get the coordinate bias
        *  bias X     bias Y      bias Z     bias N      bias E      bias U/ht
        *  0.92873	-1.83065	-0.83925	-0.63115	-1.23223	-1.73213    mean
        *  0.03027	 0.04223	 0.01276	 0.02500	 0.02875	 0.05076    std
        *  transform to WGS84(G730)(2024.0), then get the diff
stn, bias X,    bias Y ,   bias Z  ,    bias N,     bias E,     bias U       
KALU,0.9422,	-1.8691,	-0.7446,	-0.5475,	-1.2536,	-1.7506
CLMB,0.9285,	-1.8215,	-0.7179,	-0.5171,	-1.2343,	-1.7043
KEGA,0.8626,	-1.9400,	-0.7143,	-0.4853,	-1.1756,	-1.8440
MADA,0.9073,	-1.8735,	-0.7151,	-0.4896,	-1.2238,	-1.7627
mean,0.9102,	-1.8760,	-0.7230,	-0.5099,	-1.2218,	-1.7654
std ,0.0301,	 0.0422,	 0.0126,	 0.0249,	 0.0287,	 0.0503
        */
#if 0
        /* output ITRF1991(1994.0) */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 1994;
        sprintf(buffer, "WGS84(G730)(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to WGS84(G730) */
        convert_itrf2020_to_itrf1991(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
#else
        /* output WGS84 based on fitting */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 2024.0;
        sprintf(buffer, "WGS84(G730)(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to WGS84(G730)(2024.0) */
        convert_itrf2020_to_itrf1991(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
        xyz_regional[0] += 0.9102;
        xyz_regional[1] +=-1.8760;
        xyz_regional[2] +=-0.7230;
#endif
    }
    else if (code.find("TWN") != std::string::npos) /* ITRF2020(2025.0) */
    {
        /* output ITRF2020(2025.0) => change to 2025.0 (https://egnss.nlsc.gov.tw/HotNews.aspx) https://egnss.nlsc.gov.tw/content.aspx?i=20150625102013300 */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 2025.0;
        sprintf(buffer, "ITRF2020(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* predict to ITRF2020 2023.0 epoch */
        xyz_regional[0] = xyz_itrf2020[0] + vxyz_itrf2020[0] * (epoch_regional - epoch_itrf2020);
        xyz_regional[1] = xyz_itrf2020[1] + vxyz_itrf2020[1] * (epoch_regional - epoch_itrf2020);
        xyz_regional[2] = xyz_itrf2020[2] + vxyz_itrf2020[2] * (epoch_regional - epoch_itrf2020);
        vxyz_regional[0] = vxyz_itrf2020[0];
        vxyz_regional[1] = vxyz_itrf2020[1];
        vxyz_regional[2] = vxyz_itrf2020[2];
    }
    else if (code.find("THA") != std::string::npos) /* ITRF2014(2010) */
    {
        /* output ITRF2014(2010.0) */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 2010.0;
        sprintf(buffer, "ITRF2014(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2014 */
        convert_itrf2020_to_itrf2014(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("KOR") != std::string::npos) /* KGD2002=ITRF2000(2002.0) */
    {   /* https://unstats.un.org/unsd/geoinfo/rcc/docs/rccap19/ip/E_Conf.102_IP17_Korea_19th_UNRCC-AP_Session3_final.pdf */
        /* output KGD2002=ITRF2000(2002.0) */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 2000.0;
        sprintf(buffer, "KGD2002(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2000 */
        convert_itrf2020_to_itrf2000(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (code.find("MYS") != std::string::npos) /* MGRF2000=ITRF2000(2000.0), MGRF2020=ITRF2020(2020.0) */
    {
#if 0        
        /* output MGRF2000=ITRF2000(2000.0) */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 2000.0;
        sprintf(buffer, "MGRF2000(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2000 */
        convert_itrf2020_to_itrf2000(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
#else
        /* output MGRF2020=ITRF2020(2020.0) */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 2020.0;
        sprintf(buffer, "MGRF2020(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* predict to ITRF2020 2020.0 epoch */
        xyz_regional[0] = xyz_itrf2020[0] + vxyz_itrf2020[0] * (epoch_regional - epoch_itrf2020);
        xyz_regional[1] = xyz_itrf2020[1] + vxyz_itrf2020[1] * (epoch_regional - epoch_itrf2020);
        xyz_regional[2] = xyz_itrf2020[2] + vxyz_itrf2020[2] * (epoch_regional - epoch_itrf2020);
        vxyz_regional[0] = vxyz_itrf2020[0];
        vxyz_regional[1] = vxyz_itrf2020[1];
        vxyz_regional[2] = vxyz_itrf2020[2];
#endif        
    }    
    else if (code.find("ARE") != std::string::npos) /* MTRF2000=ITRF2000(2004.0) */
    {
        /* output MTRF2000=ITRF2000(2004.0) */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 2004.0;
        sprintf(buffer, "MTRF2000(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2000 */
        convert_itrf2020_to_itrf2000(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }     
    else if (is_sa_station(code)) /* SIRGAS2000=ITRF2000(2000.4) */
    {
        /* output SIRGAS2000=ITRF2000(2000.4) */
        epoch_regional = !vel_flag ? epoch_itrf2020 : 2000.4;
        sprintf(buffer, "SIRGAS2000(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ITRF2000 */
        convert_itrf2020_to_itrf2000(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (is_na_station(code))
    {
        /* output NAD83(2011)(2010.0) */
        epoch_regional = vel_flag ? 2010.0 : epoch_itrf2020;
        sprintf(buffer, "NAD83(2011)(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to NAD83(2011)(2010.0) */
        convert_itrf2020_to_nad_2011(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else if (is_eu1_station(code) || is_eu2_station(code))   /* ETRF2000(2010.0) */
    {
        /* output ETRF2000(2010.0) */
        epoch_regional = vel_flag ? 2010.0 : epoch_itrf2020;
        sprintf(buffer, "ETRF2000(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to ETRF2000 */
        convert_itrf2020_to_etrf2000(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }
    else  /* WGS84(G2139)(2025.50) */
    {
        /* output WGS84(G2139)(2025.50) */
        epoch_regional = vel_flag ? floor(epoch_itrf2020) + 0.5 : epoch_itrf2020;
        sprintf(buffer, "WGS84(G2139)(%7.2f)", epoch_regional);
        coord_name_regional = std::string(buffer);
        /* convert ITRF2020 to WGS84(G2139) */
        convert_itrf2020_to_itrf2014(xyz_itrf2020, vxyz_itrf2020, epoch_itrf2020, epoch_regional, xyz_regional, vxyz_regional);
    }

    printf("SOLU ,%s,%14s,%s,%s,%17s,%10.5f,%14.4f,%14.4f,%14.4f,%6.2f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%s,%20s,%10.5f,%14.4f,%14.4f,%14.4f,%10.6f,%10.6f,%10.6f,%17s,%10.5f,%14.4f,%14.4f,%14.4f,%10.6f,%10.6f,%10.6f\n"
        , code.c_str(), name.c_str(), stime.c_str(), ctime.c_str()
        , coord_name_itrf2020.c_str(), epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], amb_fix_rate, sigma95_xyz[0], sigma95_xyz[1], sigma95_xyz[2]
        , vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2], fname.c_str()
        , coord_name_regional.c_str(), epoch_regional, xyz_regional[0], xyz_regional[1], xyz_regional[2], vxyz_regional[0], vxyz_regional[1], vxyz_regional[2], coord_name_wgs84.c_str(), epoch_wgs84, xyz_wgs84[0], xyz_wgs84[1], xyz_wgs84[2], vxyz_wgs84[0], vxyz_wgs84[1], vxyz_wgs84[2]);

    ret = 1;

    return ret;
}
int coord_t::read_from_opus_file(const char* opusfname)
{
    int ret = 0;
    if (!opusfname) return ret;

    char buffer[255] = { 0 };

    ctime = get_current_date_string();

    FILE* fLOG = fopen(opusfname, "rt");

    type = "OPUS";

    int numofline = 0;
    int yy = 0, mm = 0, dd = 0, hour = 0, min = 0, sec = 0;
    char* val[MAXFIELD];
    const char* tempfname = strrchr(opusfname, '\\');
    if (tempfname)
    {
        strcpy(buffer, tempfname + 1);
    }
    else if (tempfname = strrchr(opusfname, '/'))
    {
        strcpy(buffer, tempfname + 1);
    }
    else
    {
        strcpy(buffer, opusfname);
    }
    char* temp = strrchr(buffer, '.');
    if (temp) temp[0] = '\0';

    fname = std::string(buffer);

    int num = parse_fields(buffer, val, '_');
    if (num >= 6)
    {
        name = std::string(val[num - 6]);
    }
    else
    {
        name = fname;
    }
    while (fLOG && !feof(fLOG) && fgets(buffer, sizeof(buffer), fLOG))
    {
        //printf("%s",buffer);
        if (temp = strchr(buffer, '\n')) temp[0] = '\0';
        if (temp = strchr(buffer, '\r')) temp[0] = '\0';
        if ((temp = strstr(buffer, "FILE:")) && !strstr(buffer, "RINEX FILE:") && !strstr(buffer, "NAV FILE:"))
        {
            remove_lead(temp + 5);
            char* temp11 = strchr(temp + 5, ' ');
            if (temp11) temp11[0] = '\0';
            temp11 = strchr(temp + 5, '.');
            if (temp11) temp11[0] = '\0';
            num = parse_fields(temp + 5, val, '_');
            if (num > 1)
            {
                name = std::string(val[1]);
            }
            else
            {
                name = std::string(val[0]);
            }
            //fname = std::string(temp + 5);
        }
        else if ((temp = strstr(buffer, "DATE:")))
        {
            remove_lead(temp + 5);
            char* temp1 = strchr(temp + 5, ',');
            if (temp1)
            {
                yy = atoi(temp1 + 1);
                char* temp2 = strchr(temp + 5, ' ');
                if (temp2)
                {
                    dd = atoi(temp2 + 1);
                    mm = 0;
                    if (strstr(temp + 5, "Jan"))
                    {
                        mm = 1;
                    }
                    else if (strstr(temp + 5, "Feb"))
                    {
                        mm = 2;
                    }
                    else if (strstr(temp + 5, "Mar"))
                    {
                        mm = 3;
                    }
                    else if (strstr(temp + 5, "Apr"))
                    {
                        mm = 4;
                    }
                    else if (strstr(temp + 5, "May"))
                    {
                        mm = 5;
                    }
                    else if (strstr(temp + 5, "Jun"))
                    {
                        mm = 6;
                    }
                    else if (strstr(temp + 5, "Jul"))
                    {
                        mm = 7;
                    }
                    else if (strstr(temp + 5, "Aug"))
                    {
                        mm = 8;
                    }
                    else if (strstr(temp + 5, "Sep"))
                    {
                        mm = 9;
                    }
                    else if (strstr(temp + 5, "Oct"))
                    {
                        mm = 10;
                    }
                    else if (strstr(temp + 5, "Nov"))
                    {
                        mm = 11;
                    }
                    else if (strstr(temp + 5, "Dec"))
                    {
                        mm = 12;
                    }
                    else
                    {
                        printf("wrong compute date %s %s\n", name.c_str(), opusfname);
                    }
                    if (mm > 0)
                    {
                        char ctime_s[20] = { 0 };
                        sprintf(ctime_s, "%04i-%02i-%02i", yy, mm, dd);
                        ctime = std::string(ctime_s);
                    }
                }
            }
        }
        else if ((temp = strstr(buffer, "START:")))
        {
            remove_lead(temp + 6);
            num = parse_fields(temp + 6, val, '/');
            if (num > 2)
            {
                char stime_s[20] = { 0 };
                sprintf(stime_s, "%04i-%02i-%02i", atoi(val[0]), atoi(val[1]), atoi(val[2]));
                stime = std::string(stime_s);
            }
        }
        else if ((temp = strstr(buffer, "STOP:")))
        {
            remove_lead(temp + 5);
            num = parse_fields(temp + 5, val, '/');
            if (num > 2)
            {
                char etime_s[20] = { 0 };
                sprintf(etime_s, "%04i-%02i-%02i", atoi(val[0]), atoi(val[1]), atoi(val[2]));
                etime = std::string(etime_s);
            }
        }
        else if ((temp = strstr(buffer, "# FIXED AMB:")))
        {
            num = sscanf(temp + 30, "%lf", &amb_fix_rate);
        }
        else if ((temp = strstr(buffer, "OVERALL RMS:")))
        {
            num = sscanf(temp + 12, "%lf", &rms);
        }
        else if (strstr(buffer, "REF FRAME:"))
        {
            temp = strstr(buffer + 55, "(EPOCH:");
            if (temp)
            {
                epoch_itrf2014 = epoch_itrf2020 = atof(temp + 7);
                temp[0] = '\0';
            }
            remove_lead(buffer + 55);
            if (strstr(buffer + 55, "ITRF2014"))
            {
                coord_name_itrf2014 = std::string(buffer + 55); buffer[55] = '\0';
            }
            else
            {
                coord_name_itrf2020 = std::string(buffer + 55); buffer[55] = '\0';
            }
            temp = strstr(buffer + 11, "(EPOCH:");
            if (temp)
            {
                epoch_regional = atof(temp + 7);
                temp[0] = '\0';
            }
            remove_lead(buffer + 11);
            if (strstr(buffer + 11, "NAD_83(2011)"))
            {
                coord_name_regional = "NAD83(2011)";
            }
            else if (strstr(buffer + 11, "NAD_83(PA11)"))
            {
                coord_name_regional = "NAD83(PA11)";
            }
            else if (strstr(buffer + 11, "NAD_83(MA11)"))
            {
                coord_name_regional = "NAD83(MA11)";
            }
            else
            {
                coord_name_regional = std::string(buffer + 11);
            }
            sprintf(buffer, "(%7.2f)", epoch_itrf2020);
            if (coord_name_itrf2014.length()>0)
                coord_name_itrf2014 += std::string(buffer);
            else 
                coord_name_itrf2020 += std::string(buffer);
            sprintf(buffer, "(%7.2f)", epoch_regional);
            coord_name_regional += std::string(buffer);
        }
        else if ((temp = strstr(buffer, "X:")) && strstr(buffer, "(m)"))
        {
            xyz_regional[0] = atof(temp + 2);
            if (coord_name_itrf2014.length()>0)
                xyz_itrf2014[0] = atof(temp + 33);
            else
                xyz_itrf2020[0] = atof(temp + 33);
            sigma95_xyz[0] = 2.0 * atof(temp + 22);
        }
        else if ((temp = strstr(buffer, "Y:")) && strstr(buffer, "(m)"))
        {
            xyz_regional[1] = atof(temp + 2);
            if (coord_name_itrf2014.length() > 0)
                xyz_itrf2014[1] = atof(temp + 33);
            else
                xyz_itrf2020[1] = atof(temp + 33);
            sigma95_xyz[1] = 2.0 * atof(temp + 22);
        }
        else if ((temp = strstr(buffer, "Z:")) && strstr(buffer, "(m)"))
        {
            xyz_regional[2] = atof(temp + 2);
            if (coord_name_itrf2014.length() > 0)
                xyz_itrf2014[2] = atof(temp + 33);
            else
                xyz_itrf2020[2] = atof(temp + 33);
            sigma95_xyz[2] = 2.0 * atof(temp + 22);
        }
        else if ((temp = strstr(buffer, "LAT:")) && strstr(buffer, "(m)"))
        {
            blh_regional[0] = (atof(temp + 4) + atof(temp + 9) / 60.0 + atof(temp + 12) / 3600.0) * D2R;
            if (coord_name_itrf2014.length() > 0)
                blh_itrf2014[0] = (atof(temp + 35) + atof(temp + 45) / 60.0 + atof(temp + 48) / 3600.0) * D2R;
            else
                blh_itrf2020[0] = (atof(temp + 35) + atof(temp + 45) / 60.0 + atof(temp + 48) / 3600.0) * D2R;
            sigma95_NEU[0] = 2.0 * atof(temp + 24);
        }
        else if ((temp = strstr(buffer, "LON:")) && strstr(buffer, "(m)") && strstr(buffer, "E"))
        {
            blh_regional[1] = (atof(temp + 4) + atof(temp + 9) / 60.0 + atof(temp + 12) / 3600.0) * D2R;
            if (coord_name_itrf2014.length() > 0)
                blh_itrf2014[1] = (atof(temp + 35) + atof(temp + 45) / 60.0 + atof(temp + 48) / 3600.0) * D2R;
            else
                blh_itrf2020[1] = (atof(temp + 35) + atof(temp + 45) / 60.0 + atof(temp + 48) / 3600.0) * D2R;
            sigma95_NEU[1] = 2.0 * atof(temp + 24);
        }
        else if ((temp = strstr(buffer, "HGT:")) && strstr(buffer, "(m)") && strstr(buffer, "EL"))
        {
            blh_regional[2] = atof(temp + 4);
            if (coord_name_itrf2014.length() > 0)
                blh_itrf2014[2] = atof(temp + 35);
            else
                blh_itrf2020[2] = atof(temp + 35);
            sigma95_NEU[2] = 2.0 * atof(temp + 24);
        }
    }

    int is_itrf2014 = !(fabs(xyz_itrf2014[0]) < 0.1 || fabs(xyz_itrf2014[1]) < 0.1 || fabs(xyz_itrf2014[2]) < 0.1 || fabs(blh_itrf2014[0]) < 0.00001 || fabs(blh_itrf2014[1]) < 0.000001 || epoch_itrf2020 < 0.001);
    int is_itrf2020 = !(fabs(xyz_itrf2020[0]) < 0.1 || fabs(xyz_itrf2020[1]) < 0.1 || fabs(xyz_itrf2020[2]) < 0.1 || fabs(blh_itrf2020[0]) < 0.00001 || fabs(blh_itrf2020[1]) < 0.000001 || epoch_itrf2020 < 0.001);
    int is_regional = !(fabs(xyz_regional[0]) < 0.1 || fabs(xyz_regional[1]) < 0.1 || fabs(xyz_regional[2]) < 0.1 || epoch_regional < 0.001);

    if (is_itrf2014 && !is_itrf2020)
    {
        double vxyz[3] = { 0 }, xyz2[3] = { 0 }, vxyz2[3] = { 0 };
        convert_itrf2014_to_itrf2020(xyz_itrf2014, vxyz, epoch_itrf2020, epoch_itrf2020, xyz_itrf2020, vxyz2);
        sprintf(buffer, "ITRF2020(%7.2f)", epoch_itrf2020);
        coord_name_itrf2020 = std::string(buffer);
        is_itrf2020 = 1;
    }
    else if (!is_itrf2014 && is_itrf2020)
    {
        double vxyz[3] = { 0 }, xyz2[3] = { 0 }, vxyz2[3] = { 0 };
        convert_itrf2020_to_itrf2014(xyz_itrf2020, vxyz, epoch_itrf2020, epoch_itrf2020, xyz_itrf2014, vxyz2);
        sprintf(buffer, "ITRF2014(%7.2f)", epoch_itrf2020);
        coord_name_itrf2014 = std::string(buffer);
        is_itrf2014 = 1;
    }
    if (is_itrf2020)
    {

        //ecef2pos_(xyz_itrf2020, blh_itrf2020, RE_GRS80, FE_GRS80);

        //get_station_code(name.c_str(), blh_itrf2020[0] * R2D, blh_itrf2020[1] * R2D, code);

        //if (get_station_vel(name.c_str(), blh_itrf2020[0] * R2D, blh_itrf2020[1] * R2D, vxyz_itrf2020))
        {
        }

        double rms3D = sqrt(sigma95_xyz[0] * sigma95_xyz[0] + sigma95_xyz[1] * sigma95_xyz[1] + sigma95_xyz[2] * sigma95_xyz[2]);

#if 0
        printf("OPUS ,%s,%14s,%s,%s,%17s,%10.5f,%14.4f,%14.4f,%14.4f,%6.2f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%s,%20s,%10.5f,%14.4f,%14.4f,%14.4f,%17s,%10.5f,%14.4f,%14.4f,%14.4f,%10.6f\n"
            , code.c_str(), name.c_str(), stime.c_str(), ctime.c_str()
            , coord_name_itrf2020.c_str(), epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], amb_fix_rate, sigma95_xyz[0], sigma95_xyz[1], sigma95_xyz[2]
            , vxyz_itrf2020[0], vxyz_itrf2020[1], vxyz_itrf2020[2], fname.c_str()
            , coord_name_regional.c_str(), epoch_regional, xyz_regional[0], xyz_regional[1], xyz_regional[2], coord_name_itrf2014.c_str(), epoch_itrf2020, xyz_itrf2014[0], xyz_itrf2014[1], xyz_itrf2014[2], rms);
#endif
        printf("OPUS ,%14s,%s,%s,%17s,%10.5f,%14.4f,%14.4f,%14.4f,%6.2f,%10.6f,%10.6f,%10.6f,%20s,%10.5f,%14.4f,%14.4f,%14.4f,%10.4f,%10.4f,%s\n"
            , name.c_str(), stime.c_str(), ctime.c_str()
            , coord_name_itrf2020.c_str(), epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], amb_fix_rate, sigma95_xyz[0], sigma95_xyz[1], sigma95_xyz[2]
            , coord_name_regional.c_str(), epoch_regional, xyz_regional[0], xyz_regional[1], xyz_regional[2], rms, rms3D, fname.c_str());

        ret = 1;

    }
    if (fLOG) fclose(fLOG);
    return ret;
}
int coord_t::read_from_nrcan_file(const char* nrcanfname)
{
    int ret = 0;
    if (!nrcanfname) return ret;
    FILE* fLOG = fopen(nrcanfname, "rt");
    char buffer[255] = { 0 };

    ctime = get_current_date_string();

    int numofline = 0;
    int yy = 0, mm = 0, dd = 0, hour = 0, min = 0, sec = 0, doy = 0;
    char* temp = strrchr(buffer, '-');

    char* val[MAXFIELD];

    double blh_itrf2020[3] = { 0 };

    const char* temp3 = strrchr(nrcanfname, '\\');
    if (temp3)
    {
        strcpy(buffer, temp3 + 1);
    }
    else if (temp3 = strrchr(nrcanfname, '/'))
    {
        strcpy(buffer, temp3 + 1);
    }
    else
    {
        strcpy(buffer, nrcanfname);
    }
    temp = strrchr(buffer, '.');
    if (temp) temp[0] = '\0';

    fname = std::string(buffer);

    int num = parse_fields(buffer, val, '_');
    if (num >= 5)
    {
        name = std::string(val[num - 5]);
    }
    else
    {
        name = fname;
    }

    type = "NRCAN";

    while (fLOG && !feof(fLOG))
    {
        fgets(buffer, sizeof(buffer), fLOG);
        temp = strrchr(buffer, '\n');
        if (temp)
            temp[0] = '\0';
        temp = strrchr(buffer, '\r');
        if (temp)
            temp[0] = '\0';
        if (strlen(buffer) >= 4 && buffer[0] == 'R' && buffer[1] == 'N' && buffer[2] == 'X')
        {
            temp = strrchr(buffer + 4, '.');
            if (temp) temp[0] = '\0';
            //fname = std::string(buffer + 4);
        }
        if (strlen(buffer) >= 4 && buffer[0] == 'N' && buffer[1] == 'O' && buffer[2] == 'W')
        {
            ctime = std::string(buffer + 4);
            std::size_t nloc = ctime.find(' ');
            if (nloc != std::string::npos)
            {
                ctime = ctime.substr(0, nloc);
            }
        }
        if (strlen(buffer) >= 4 && buffer[0] == 'B' && buffer[1] == 'E' && buffer[2] == 'G')
        {
            stime = std::string(buffer + 4);
            std::size_t nloc = stime.find(' ');
            if (nloc != std::string::npos)
            {
                stime = stime.substr(0, nloc);
            }
        }
        if (strlen(buffer) >= 4 && buffer[0] == 'E' && buffer[1] == 'N' && buffer[2] == 'D')
        {
            etime = std::string(buffer + 4);
            std::size_t nloc = etime.find(' ');
            if (nloc != std::string::npos)
            {
                etime = etime.substr(0, nloc);
            }
        }
        if (strlen(buffer) >= 4 && buffer[0] == 'I' && buffer[1] == 'A' && buffer[2] == 'R')
        {
            amb_fix_rate = atof(buffer + 4);
        }
        if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[6] == 'X')
        {
            coord_name_itrf2020 = std::string(buffer + 8).substr(0, 5);
            yy = atoi(buffer + 14);
            doy = atoi(buffer + 17);
            epoch_itrf2020 = 2000.0 + yy + doy / 365.0;
            xyz_itrf2020[0] = atof(buffer + 44);
            sigma95_xyz[0] = atof(buffer + 73);
        }
        if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[6] == 'Y')
        {
            coord_name_itrf2020 = std::string(buffer + 8).substr(0, 5);
            yy = atoi(buffer + 14);
            doy = atoi(buffer + 17);
            epoch_itrf2020 = 2000.0 + yy + doy / 365.0;
            xyz_itrf2020[1] = atof(buffer + 44);
            sigma95_xyz[1] = atof(buffer + 73);
        }
        if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[6] == 'Z')
        {
            coord_name_itrf2020 = std::string(buffer + 8).substr(0, 5);
            yy = atoi(buffer + 14);
            doy = atoi(buffer + 17);
            epoch_itrf2020 = 2000.0 + yy + doy / 365.0;
            xyz_itrf2020[2] = atof(buffer + 44);
            sigma95_xyz[2] = atof(buffer + 73);
        }
        if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[4] == 'L' && buffer[5] == 'A' && buffer[6] == 'T')
        {
            coord_name_itrf2020 = std::string(buffer + 8).substr(0, 5);
            yy = atoi(buffer + 14);
            doy = atoi(buffer + 17);
            epoch_itrf2020 = 2000.0 + yy + doy / 365.0;
            char dd_str[7] = { 0 };
            strncpy(dd_str, buffer + 44, 6);
            dd_str[6] = '\0';
            temp = strchr(dd_str, '-');
            int is_neg = 0;
            if (temp)
            {
                is_neg = 1;
                temp[0] = ' ';
            }
            int dd = atoi(dd_str);// atoi(buffer + 44);
            int mm = atoi(buffer + 50);
            double ss = atof(buffer + 53);
            if (is_neg)
            {
                blh_itrf2020[0] = -(dd + mm / 60.0 + ss / 3600.0);
            }
            else
            {
                blh_itrf2020[0] = (dd + mm / 60.0 + ss / 3600.0);
            }
            sigma95_NEU[0] = atof(buffer + 73);
        }
        if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[4] == 'L' && buffer[5] == 'O' && buffer[6] == 'N')
        {
            coord_name_itrf2020 = std::string(buffer + 8).substr(0, 5);
            yy = atoi(buffer + 14);
            doy = atoi(buffer + 17);
            epoch_itrf2020 = 2000.0 + yy + doy / 365.0;
            char dd_str[7] = { 0 };
            strncpy(dd_str, buffer + 44, 6);
            dd_str[6] = '\0';
            temp = strchr(dd_str, '-');
            int is_neg = 0;
            if (temp)
            {
                is_neg = 1;
                temp[0] = ' ';
            }
            int dd = atoi(dd_str);// atoi(buffer + 44);
            int mm = atoi(buffer + 50);
            double ss = atof(buffer + 53);
            if (is_neg)
            {
                blh_itrf2020[1] = -(dd + mm / 60.0 + ss / 3600.0);
            }
            else
            {
                blh_itrf2020[1] = (dd + mm / 60.0 + ss / 3600.0);
            }
            sigma95_NEU[1] = atof(buffer + 73);
        }
        if (strlen(buffer) >= 83 && buffer[0] == 'P' && buffer[1] == 'O' && buffer[2] == 'S' && buffer[4] == 'H' && buffer[5] == 'G' && buffer[6] == 'T')
        {
            coord_name_itrf2020 = std::string(buffer + 8).substr(0, 5);
            int yy = atoi(buffer + 14);
            int doy = atoi(buffer + 17);
            epoch_itrf2020 = 2000.0 + yy + doy / 365.0;
            blh_itrf2020[2] = atof(buffer + 44);
            sigma95_NEU[2] = atof(buffer + 73);
        }
    }

    if (coord_name_itrf2020.find('(') == std::string::npos)
    {
        memset(buffer, 0, sizeof(buffer));
        sprintf(buffer, "(%.2f)", epoch_itrf2020);
        coord_name_itrf2020 += std::string(buffer);
    }

    if (fabs(xyz_itrf2020[0]) < 0.001 || fabs(xyz_itrf2020[1]) < 0.001 || fabs(xyz_itrf2020[2]) < 0.001 || epoch_itrf2020 < 0.0001 || sigma95_xyz[0] < 0.0000001 || sigma95_xyz[1] < 0.0000001 || sigma95_xyz[2] < 0.0000001)
    {

    }
    else
    {
        ecef2pos_(xyz_itrf2020, blh_itrf2020, RE_GRS80, FE_GRS80);

        get_station_code(name.c_str(), blh_itrf2020[0] * R2D, blh_itrf2020[1] * R2D, code);

        sprintf(buffer, "ITRF2020(%7.2f)", epoch_itrf2020);
        coord_name_itrf2020 = std::string(buffer);

        printf("NRCAN,%s,%14s,%s,%s,ITRF2020(%7.2f),%10.5f,%14.4f,%14.4f,%14.4f,%6.2f,%10.6f,%10.6f,%10.6f,%s\n"
            , code.c_str(), name.c_str(), stime.c_str(), ctime.c_str(), epoch_itrf2020, epoch_itrf2020, xyz_itrf2020[0], xyz_itrf2020[1], xyz_itrf2020[2], amb_fix_rate, sigma95_xyz[0], sigma95_xyz[1], sigma95_xyz[2], fname.c_str());

        ret = 1;

    }
    if (fLOG) fclose(fLOG);
    return ret;
}

int coord_t::read_estimation_list(char* buffer)
{
    /* export from estimation list
    Station,Date,Platform,Status,Result,Message,Rejected Epochs,Fixed Ambiguities,Sigmas(95%) North,SigmasE(95%) East,SigmasU(95%) Height,Primitive,ITRF2020,WGS84,Regional
    */
    int ret = 0;
    if (strstr(buffer, "Station") || strstr(buffer, "Platform")) return ret;
    char* temp = nullptr;
    while (temp = strstr(buffer, ":")) temp[0] = ',';
    while (temp = strstr(buffer, ";")) temp[0] = ',';
    /*
    * E465B86A0745,2025-03-16,NRCAN,Done,Success,,0,88.7,0.0029,0.0025,0.0117,name: ITRF2020(2025.21);epoch: 2025.20548;x: 1005413.6541;y: -4878972.1712;z: 3970627.8834,name: ITRF2020(2015.00);epoch: 2015;x: 1005413.8048;y: -4878972.1658;z: 3970627.8574,name: WGS84(G2139)(2025.50);epoch: 2025.5;x: 1005413.6479;y: -4878972.1713;z: 3970627.886,name: NAD83(2011)(2010.00);epoch: 2010;x: 1005414.632;y: -4878973.6103;z: 3970627.9347
    */

    char* val[100];
    int num = parse_fields(buffer, val, ',');
    if (num > 50 && strstr(val[4], "Success"))
    {
        name = std::string(val[0]);
        stime = std::string(val[1]);
        type = std::string(val[2]);
        amb_fix_rate = atof(val[7]);
        sigma95_NEU[0] = atof(val[8]);
        sigma95_NEU[1] = atof(val[9]);
        sigma95_NEU[2] = atof(val[10]);
        coord_name_itrf2020 = std::string(val[12]);
        epoch_itrf2020 = atof(val[14]);
        xyz_itrf2020[0] = atof(val[16]);
        xyz_itrf2020[1] = atof(val[18]);
        xyz_itrf2020[2] = atof(val[20]);
        coord_name_itrf2020_2015 = std::string(val[22]);
        epoch_itrf2020_2015 = atof(val[24]);
        xyz_itrf2020_2015[0] = atof(val[26]);
        xyz_itrf2020_2015[1] = atof(val[28]);
        xyz_itrf2020_2015[2] = atof(val[30]);
        coord_name_wgs84 = std::string(val[32]);
        epoch_wgs84 = atof(val[34]);
        xyz_wgs84[0] = atof(val[36]);
        xyz_wgs84[1] = atof(val[38]);
        xyz_wgs84[2] = atof(val[40]);
        coord_name_regional = std::string(val[42]);
        epoch_regional = atof(val[44]);
        xyz_regional[0] = atof(val[46]);
        xyz_regional[1] = atof(val[48]);
        xyz_regional[2] = atof(val[50]);
        ret = 1;
    }
    return ret;
}

int coord_t::read_history_list(char* buffer)
{
    int ret = 0;
    if (strstr(buffer, "Station") || strstr(buffer, "Platform")) return ret;
    char* temp = nullptr;
    //while (temp = strstr(buffer, ":")) temp[0] = ',';
    //while (temp = strstr(buffer, ";")) temp[0] = ',';
    /*
    * Station,Date(GNSS Time),Computing time(UTC Time),Platform,Rejected Epochs,Fixed Ambiguities,Sigmas(95%) North,Sigmas(95%) East,Sigmas(95%) Height,ITRF2020 name,epoch,x,y,z,sigma95X,sigma95Y,sigma95Z,ITRF2015 name,epoch,x,y,z,vx,vy,vz,WGS84 name,x,y,z,vx,vy,vz,Regional name,epoch,x,y,z,vx,vy,vz
    * G001,2025-03-15,2025-03-16 01:02:35,NRCAN,0,99.09,0.0028,0.0023,0.0096,ITRF2020(2025.20),2025.20274,-2687303.0419,-4302926.7293,3852731.27,0.0047,0.007,0.0059,ITRF2020(2015.00),2015,-2687302.8053,-4302926.9705,3852731.1716,-0.023185,0.023645,0.009641,WGS84(G2139)(2025.50),2025.5,-2687303.0491,-4302926.7224,3852731.2747,-0.023185,0.023545,0.009841,NAD83(2011)(2010.00),2010,-2687301.9056,-4302928.3555,3852731.1462,-0.006983,0.023831,0.019383
    * G001,2025-03-01,2025-03-02 00:57:32,NRCAN,0,97.92,0.0028,0.0023,0.0096,ITRF2020(2025.2),2025.16,-2687303.041,-4302926.7311,3852731.2662,0.0047,0.007,0.0059,ITRF2020(2015.00),2015,-2687302.805,-4302926.971,3852731.168,-23.19,23.65,9.64,WGS84(G2139)(2025.50),2025.5,-2687303.049,-4302926.723,3852731.271,0,0,0,NAD83(2011)(2010.00),2010,-2687301.906,-4302928.356,3852731.143,0,0,0
    */
    char* val[100];
    int num = parse_fields(buffer, val, ',');
    if (num > 40)
    {
        name = std::string(val[0]);
        stime = std::string(val[1]);
        ctime = std::string(val[2]);
        type = std::string(val[3]);
        amb_fix_rate = atof(val[5]);
        sigma95_NEU[0] = atof(val[6]);
        sigma95_NEU[1] = atof(val[7]);
        sigma95_NEU[2] = atof(val[8]);
        coord_name_itrf2020 = std::string(val[9]);
        epoch_itrf2020 = atof(val[10]);
        xyz_itrf2020[0] = atof(val[11]);
        xyz_itrf2020[1] = atof(val[12]);
        xyz_itrf2020[2] = atof(val[13]);
        sigma95_xyz[0] = atof(val[14]);
        sigma95_xyz[1] = atof(val[15]);
        sigma95_xyz[2] = atof(val[16]);
        coord_name_itrf2020_2015 = std::string(val[17]);
        epoch_itrf2020_2015 = atof(val[18]);
        xyz_itrf2020_2015[0] = atof(val[19]);
        xyz_itrf2020_2015[1] = atof(val[20]);
        xyz_itrf2020_2015[2] = atof(val[21]);
        vxyz_itrf2020[0] = atof(val[22]);
        vxyz_itrf2020[1] = atof(val[23]);
        vxyz_itrf2020[2] = atof(val[24]);
        double vxyz3d = sqrt(vxyz_itrf2020[0] * vxyz_itrf2020[0] + vxyz_itrf2020[1] * vxyz_itrf2020[1] + vxyz_itrf2020[2] * vxyz_itrf2020[2]);
        if (vxyz3d > 0.5) /* mm/year */
        {
            vxyz_itrf2020[0] /= 1000.0;
            vxyz_itrf2020[1] /= 1000.0;
            vxyz_itrf2020[2] /= 1000.0;
        }
        coord_name_wgs84 = std::string(val[25]);
        epoch_wgs84 = atof(val[26]);
        xyz_wgs84[0] = atof(val[27]);
        xyz_wgs84[1] = atof(val[28]);
        xyz_wgs84[2] = atof(val[29]);
        vxyz_wgs84[0] = atof(val[30]);
        vxyz_wgs84[1] = atof(val[31]);
        vxyz_wgs84[2] = atof(val[32]);
        if (vxyz3d > 0.5) /* mm/year */
        {
            vxyz_wgs84[0] /= 1000.0;
            vxyz_wgs84[1] /= 1000.0;
            vxyz_wgs84[2] /= 1000.0;
        }
        coord_name_regional = std::string(val[33]);
        epoch_regional = atof(val[34]);
        xyz_regional[0] = atof(val[35]);
        xyz_regional[1] = atof(val[36]);
        xyz_regional[2] = atof(val[37]);
        vxyz_regional[0] = atof(val[38]);
        vxyz_regional[1] = atof(val[39]);
        vxyz_regional[2] = atof(val[40]);
        if (vxyz3d > 0.5) /* mm/year */
        {
            vxyz_regional[0] /= 1000.0;
            vxyz_regional[1] /= 1000.0;
            vxyz_regional[2] /= 1000.0;
        }
        ret = 1;
    }
    return ret;
}