
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <map>
#include <algorithm>
#include <filesystem>

#include "coord.h"
#include "coord_sol.h"

#ifdef _PLOT_
/* need to install matplotlib-cpp at upper directory https://github.com/yydgis/matplotlib-cpp.git */
/* need to install python and add the path => normally go to cmd, and use path command to locate the python path, for example C:\Users\yydgi\AppData\Local\Programs\Python\Python310 */
/* need to install numpy and add the path for example, C:\Users\yydgi\AppData\Local\Programs\Python\Python310\Lib\site-packages\numpy\_core\include */
/* add python312.lib in the path and dependencies, note this only works for release build with x64 platform, need to find the other lib for debug version and win32 version */

#include "../matplotlib-cpp/matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;
#endif

/* help text -----------------------------------------------------------------*/
static const char* help[] = {
" coordinate ",
" coordinate update_opus_file opus_sol_file      => update OPUS solution",
" coordinate update_opus_list opus_sol_list      => update OPUS solution in a batch",
" coordinate update_nrcan_file nrcan_sol_file    => update NRCAN solution",
" coordinate update_nrcan_list nrcan_sol_list    => update NRCAN solution in a batch",
" coordinate update_solu_file fname              => update solution in a file",
" coordinate check_station_solu fname            => check the station solutions from solution database",
" coordinate check_estimate_list fname           => check the station estimation list",
};

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i = 0; i < (int)(sizeof(help) / sizeof(*help)); i++) 
    {
        if (i==0)
            fprintf(stderr, "%2c %s\n", ' ', help[i]);
        else
            fprintf(stderr, "%02i %s\n", i, help[i]);
    }
    exit(0);
}

// Function to trim leading whitespace
void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
        }));
}

// Function to trim trailing whitespace
void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
        }).base(), s.end());
}

// Function to trim both leading and trailing whitespace
void trim(std::string& s) {
    ltrim(s);
    rtrim(s);
}

static int read_ccsv_file_all(const char* fname, std::map<std::string, std::vector<coord_t>>& mCoords)
{
    int ret = 0;
    FILE* fCCSV = fopen(fname, "r");
    char buffer[512] = { 0 };
    char* temp = nullptr;
    char* val[100];

    while (fCCSV && !feof(fCCSV) && fgets(buffer, sizeof(buffer), fCCSV))
    {
        /* Station,Date,Platform,Status,Result,Message,Rejected Epochs,Fixed Ambiguities,Sigmas(95%) North,SigmasE(95%) East,SigmasU(95%) Height,Primitive,ITRF2020,WGS84,Regional       */
        if (strstr(buffer, "Station") && strstr(buffer, "Date")) continue;
        if (temp = strchr(buffer, '\n')) temp[0] = '\0';
        int num = parse_fields(buffer, val, ',');
        if (num < 41) continue;
        coord_t coord;
        double blh[3] = { 0 };
        coord.name = std::string(val[0]);
        coord.stime = std::string(val[1]);
        coord.amb_fix_rate = atof(val[5]);
        double sigmaN = atof(val[6]);
        double sigmaE = atof(val[7]);
        double sigmaU = atof(val[8]);
        coord.sigma95_NEU[0] = sigmaN;
        coord.sigma95_NEU[1] = sigmaE;
        coord.sigma95_NEU[2] = sigmaU;
        coord.coord_name_itrf2020 = std::string(val[9]);
        coord.epoch_itrf2020 = atof(val[10]);
        coord.xyz_itrf2020[0] = atof(val[11]);
        coord.xyz_itrf2020[1] = atof(val[12]);
        coord.xyz_itrf2020[2] = atof(val[13]);
        coord.vxyz_itrf2020[0] = atof(val[22]);
        coord.vxyz_itrf2020[1] = atof(val[23]);
        coord.vxyz_itrf2020[2] = atof(val[24]);
        ecef2pos_(coord.xyz_itrf2020, blh, RE_GRS80, FE_GRS80);

        trim(coord.coord_name_itrf2020);
        trim(coord.name);

        double C_en[3][3] = { 0 };
        double lat = blh[0];
        double lon = blh[1];

        C_en[0][0] = -sin(lat) * cos(lon);
        C_en[1][0] = -sin(lat) * sin(lon);
        C_en[2][0] = cos(lat);
        C_en[0][1] = -sin(lon);
        C_en[1][1] = cos(lon);
        C_en[2][1] = 0.0;
        C_en[0][2] = -cos(lat) * cos(lon);
        C_en[1][2] = -cos(lat) * sin(lon);
        C_en[2][2] = -sin(lat);

        /* dXYZ = C_en*dNED */

        /* cov(xyz) = C_en*cov(ned)*C_en' */
        double covX = C_en[0][0] * sigmaN * sigmaN * C_en[0][0] + C_en[0][1] * sigmaE * sigmaE * C_en[0][1] + C_en[0][2] * sigmaU * sigmaU * C_en[0][2];
        double covY = C_en[1][0] * sigmaN * sigmaN * C_en[1][0] + C_en[1][1] * sigmaE * sigmaE * C_en[1][1] + C_en[1][2] * sigmaU * sigmaU * C_en[1][2];
        double covZ = C_en[2][0] * sigmaN * sigmaN * C_en[2][0] + C_en[2][1] * sigmaE * sigmaE * C_en[2][1] + C_en[2][2] * sigmaU * sigmaU * C_en[2][2];
        coord.sigma95_xyz[0] = sqrt(covX);
        coord.sigma95_xyz[1] = sqrt(covY);
        coord.sigma95_xyz[2] = sqrt(covZ);

        std::map<std::string, std::vector<coord_t>>::iterator pCoords = mCoords.find(coord.name);
        if (pCoords == mCoords.end())
        {
            mCoords[coord.name].push_back(coord);
        }
        else
        {
            std::vector<coord_t>::iterator pCoord = std::find(pCoords->second.begin(), pCoords->second.end(), coord);
            if (pCoord == pCoords->second.end())
            {
                pCoords->second.push_back(coord);
            }
            else
            {
                *pCoord = coord;
            }
        }
    }
    if (fCCSV) fclose(fCCSV);
    return ret;
}

static int read_ccsv_file1(const char* fname, std::map<std::string, std::vector<coord_t>>& mCoords)
{
    int ret = 0;
    FILE* fCCSV = fopen(fname, "r");
    char buffer[512] = { 0 };
    char* temp = nullptr;
    char* val[100];

    while (fCCSV && !feof(fCCSV) && fgets(buffer, sizeof(buffer), fCCSV))
    {
        /* Station,Date,Platform,Status,Result,Message,Rejected Epochs,Fixed Ambiguities,Sigmas(95%) North,SigmasE(95%) East,SigmasU(95%) Height,Primitive,ITRF2020,WGS84,Regional       */
        if (strstr(buffer, "Station") && strstr(buffer, "Date")) continue;
        if (temp = strchr(buffer, '\n')) temp[0] = '\0';
        while (temp = strchr(buffer, ':')) temp[0] = ',';
        while (temp = strchr(buffer, ';')) temp[0] = ',';
        int num = parse_fields(buffer, val, ',');
        if (num < 51) continue;
        coord_t coord;
        double blh[3] = { 0 };
        coord.name = std::string(val[0]);
        coord.stime = std::string(val[1]);
        coord.amb_fix_rate = atof(val[7]);
        double sigmaN = atof(val[8]);
        double sigmaE = atof(val[9]);
        double sigmaU = atof(val[10]);
        coord.sigma95_NEU[0] = sigmaN;
        coord.sigma95_NEU[1] = sigmaE;
        coord.sigma95_NEU[2] = sigmaU;
        coord.coord_name_itrf2020 = std::string(val[12]);
        coord.epoch_itrf2020 = atof(val[14]);
        coord.xyz_itrf2020[0] = atof(val[16]);
        coord.xyz_itrf2020[1] = atof(val[18]);
        coord.xyz_itrf2020[2] = atof(val[20]);
        coord.epoch_itrf2020_2015 = atof(val[24]);
        coord.xyz_itrf2020_2015[0] = atof(val[26]);
        coord.xyz_itrf2020_2015[1] = atof(val[28]);
        coord.xyz_itrf2020_2015[2] = atof(val[30]);
        coord.vxyz_itrf2020[0] = (coord.xyz_itrf2020[0] - coord.xyz_itrf2020_2015[0]) / (coord.epoch_itrf2020 - coord.epoch_itrf2020_2015);
        coord.vxyz_itrf2020[1] = (coord.xyz_itrf2020[1] - coord.xyz_itrf2020_2015[1]) / (coord.epoch_itrf2020 - coord.epoch_itrf2020_2015);
        coord.vxyz_itrf2020[2] = (coord.xyz_itrf2020[2] - coord.xyz_itrf2020_2015[2]) / (coord.epoch_itrf2020 - coord.epoch_itrf2020_2015);
        ecef2pos_(coord.xyz_itrf2020, blh, RE_GRS80, FE_GRS80);

        trim(coord.coord_name_itrf2020);
        trim(coord.name);

        double C_en[3][3] = { 0 };
        double lat = blh[0];
        double lon = blh[1];

        C_en[0][0] = -sin(lat) * cos(lon);
        C_en[1][0] = -sin(lat) * sin(lon);
        C_en[2][0] = cos(lat);
        C_en[0][1] = -sin(lon);
        C_en[1][1] = cos(lon);
        C_en[2][1] = 0.0;
        C_en[0][2] = -cos(lat) * cos(lon);
        C_en[1][2] = -cos(lat) * sin(lon);
        C_en[2][2] = -sin(lat);

        /* dXYZ = C_en*dNED */

        /* cov(xyz) = C_en*cov(ned)*C_en' */
        double covX = C_en[0][0] * sigmaN * sigmaN * C_en[0][0] + C_en[0][1] * sigmaE * sigmaE * C_en[0][1] + C_en[0][2] * sigmaU * sigmaU * C_en[0][2];
        double covY = C_en[1][0] * sigmaN * sigmaN * C_en[1][0] + C_en[1][1] * sigmaE * sigmaE * C_en[1][1] + C_en[1][2] * sigmaU * sigmaU * C_en[1][2];
        double covZ = C_en[2][0] * sigmaN * sigmaN * C_en[2][0] + C_en[2][1] * sigmaE * sigmaE * C_en[2][1] + C_en[2][2] * sigmaU * sigmaU * C_en[2][2];
        coord.sigma95_xyz[0] = sqrt(covX);
        coord.sigma95_xyz[1] = sqrt(covY);
        coord.sigma95_xyz[2] = sqrt(covZ);
        std::map<std::string, std::vector<coord_t>>::iterator pCoords = mCoords.find(coord.name);
        if (pCoords == mCoords.end())
        {
            mCoords[coord.name].push_back(coord);
        }
        else
        {
            std::vector<coord_t>::iterator pCoord = std::find(pCoords->second.begin(), pCoords->second.end(), coord);
            if (pCoord == pCoords->second.end())
            {
                pCoords->second.push_back(coord);
            }
            else
            {
                *pCoord = coord;
            }
        }
    }
    if (fCCSV) fclose(fCCSV);
    return ret;
}

static int plot_sol_diff(std::string name, std::vector<coord_t>& mCoords, double sN, double sE, double sU)
{
#ifdef _PLOT_
    char title_buffer[255] = { 0 };
    if (name.length() > 5)
        sprintf(title_buffer, "%5s Drift in North, East and Up [m]", name.substr(name.length() - 5).c_str());
    else
        sprintf(title_buffer, "%s Drift in North, East and Up [m]", name.c_str());
    //plt::figure();
    //std::vector<std::string> dd;
    std::vector<double> dt;
    std::vector<double> dN;
    std::vector<double> dE;
    std::vector<double> dU;
    int idx = 0;
    for (std::vector<coord_t>::iterator pCoord = mCoords.begin(); pCoord != mCoords.end(); ++pCoord)
    {
        //dd.push_back(pCoord->stime);
        dt.push_back(pCoord->epoch_itrf2020 - 2000.0);
        //dt.push_back(idx++);
        dN.push_back(pCoord->dN);
        dE.push_back(pCoord->dE);
        dU.push_back(pCoord->dU);
    }

    plt::title(title_buffer);
    char name_n[255] = { 0 };
    char name_e[255] = { 0 };
    char name_u[255] = { 0 };
    sprintf(name_n, "N RMS %.4f [m]", sN);
    sprintf(name_e, "E RMS %.4f [m]", sE);
    sprintf(name_u, "U RMS %.4f [m]", sU);
    plt::named_plot(name_n, dt, dN, ".-");
    plt::named_plot(name_e, dt, dE, "+-");
    plt::named_plot(name_u, dt, dU, "*-");
    plt::legend();
    plt::grid(true);

    std::filesystem::path dir_path = "images";
    if (!std::filesystem::exists(dir_path))
    {
        std::filesystem::create_directory(dir_path);
    }
    sprintf(title_buffer, "images\\%s-drift", name.c_str());
    plt::save(title_buffer);
    //plt::show(false);
    plt::cla();
#endif
    return 0;
}

static int process_data_shift(std::map<std::string, std::vector<coord_t>>& mCoords)
{
    int ret = 0;
    printf("Station,Date,Coordinate System,amb fix rate,Sigma N[m],Sigma E[m],Sigma U[m],X [m],Y [m],Z [m],Vx[m/year],Vy[m/year],Vz[m/year],lat [deg],lon [deg], ellipsoid height [m], sN[m], sE[m], sU[m], rms3D[m], count\n");
    for (std::map<std::string, std::vector<coord_t>>::iterator pCoords = mCoords.begin(); pCoords != mCoords.end(); ++pCoords)
    {
        std::sort(pCoords->second.begin(), pCoords->second.end());
        std::vector<coord_t>::reverse_iterator pLast = pCoords->second.rbegin();
        if (fabs(pLast->xyz_itrf2020[0]) < 0.001 || fabs(pLast->xyz_itrf2020[1]) < 0.001 || fabs(pLast->xyz_itrf2020[2]) < 0.001)
        {
            continue;
        }
        double xyz0[3] = { pLast->xyz_itrf2020[0], pLast->xyz_itrf2020[1], pLast->xyz_itrf2020[2] };
        //double vxyz[3] = { pLast->vxyz_itrf2020[0], pLast->vxyz_itrf2020[1], pLast->vxyz_itrf2020[2] };
        double blh0[3] = { 0 };
        ecef2pos_(xyz0, blh0, RE_GRS80, FE_GRS80);
        double dN = 0;
        double dE = 0;
        double dU = 0;
        double xyz[3] = { 0 };
        double dxyz[3] = { 0 };
        double dned[3] = { 0 };
        int count = 0;
        double ref_epoch = pLast->epoch_itrf2020;
        double dt = 0;
        for (std::vector<coord_t>::reverse_iterator pCoord = pCoords->second.rbegin(); pCoord != pCoords->second.rend(); ++pCoord)
        {
            dt = ref_epoch - pCoord->epoch_itrf2020;
            xyz[0] = pCoord->xyz_itrf2020[0];
            xyz[1] = pCoord->xyz_itrf2020[1];
            xyz[2] = pCoord->xyz_itrf2020[2];
            //xyz[0] -= vxyz[0] * dt;
            //xyz[0] -= vxyz[1] * dt;
            //xyz[0] -= vxyz[2] * dt;
            dxyz[0] = xyz[0] - xyz0[0];
            dxyz[1] = xyz[1] - xyz0[1];
            dxyz[2] = xyz[2] - xyz0[2];
            ecef2pos_(xyz, pCoord->blh_itrf2020, RE_GRS80, FE_GRS80);
            XYZ_to_NED(blh0[0], blh0[1], dxyz, dned);
            pCoord->dN = dned[0];
            pCoord->dE = dned[1];
            pCoord->dU = -dned[2];
            dN += pCoord->dN;
            dE += pCoord->dE;
            dU += pCoord->dU;
            ++count;
        }
        if (count >= 7)
        {
            dN /= count;
            dE /= count;
            dU /= count;
            double sN = 0;
            double sE = 0;
            double sU = 0;
            for (std::vector<coord_t>::iterator pCoord = pCoords->second.begin(); pCoord != pCoords->second.end(); ++pCoord)
            {
                pCoord->dN -= dN;
                pCoord->dE -= dE;
                pCoord->dU -= dU;
                sN += pCoord->dN * pCoord->dN;
                sE += pCoord->dE * pCoord->dE;
                sU += pCoord->dU * pCoord->dU;
            }
            sN = count > 1 ? sqrt(sN / (count - 1)) : sqrt(sN / count);
            sE = count > 1 ? sqrt(sE / (count - 1)) : sqrt(sE / count);
            sU = count > 1 ? sqrt(sU / (count - 1)) : sqrt(sU / count);

            dned[0] = dN;
            dned[1] = dE;
            dned[2] = -dU;

            NED_to_XYZ(blh0[0], blh0[1], dned, dxyz);

            xyz0[0] += dxyz[0];
            xyz0[1] += dxyz[1];
            xyz0[2] += dxyz[2];

            ecef2pos_(xyz0, blh0, RE_GRS80, FE_GRS80);

            double rms3D = sqrt(sN * sN + sE * sE + sU * sU);
            printf("%12s,%s,%-17s,%7.2f,%10.6f,%10.6f,%10.6f,%14.4f,%14.4f,%14.4f,%10.6f,%10.6f,%10.6f,%14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%i\n"
                , pLast->name.c_str(), pLast->stime.c_str(), pLast->coord_name_itrf2020.c_str(), pLast->amb_fix_rate, pLast->sigma95_NEU[0], pLast->sigma95_NEU[1], pLast->sigma95_NEU[2]
                , xyz0[0], xyz0[1], xyz0[2], pLast->vxyz_itrf2020[0], pLast->vxyz_itrf2020[1], pLast->vxyz_itrf2020[2], blh0[0] * R2D, blh0[1] * R2D, blh0[2], sN, sE, sU, rms3D, count);
#if 0
            char fname[255] = { 0 };
            sprintf(fname, "%s.csv", pCoords->first.c_str());
            FILE* fOUT = fopen(fname, "w");

            if (fOUT)
                fprintf(fOUT, "Station,Date,Coordinate System,amb fix rate,Sigma N[m],Sigma E[m],Sigma U[m],X [m],Y [m],Z [m],lat [deg],lon [deg], ellipsoid height [m], dNorth[m], dEast[m], dUp[m],%14.4f,%14.4f,%14.4f,%10.6f,%10.6f,%10.6f,%14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f\n"
                    , xyz0[0], xyz0[1], xyz0[2], pLast->vxyz_itrf2020[0], pLast->vxyz_itrf2020[1], pLast->vxyz_itrf2020[2], blh0[0] * R2D, blh0[1] * R2D, blh0[2], sN, sE, sU);
            for (std::vector<coord_t>::iterator pCoord = pCoords->second.begin(); pCoord != pCoords->second.end(); ++pCoord)
            {
                if (fOUT)
                    fprintf(fOUT, "%12s,%s,%s,%7.2f,%10.6f,%10.6f,%10.6f,%14.4f,%14.4f,%14.4f,%14.9f,%14.9f,%10.4f,%10.4f,%10.4f,%10.4f\n"
                        , pCoord->name.c_str(), pCoord->stime.c_str(), pCoord->coord_name_itrf2020.c_str(), pCoord->amb_fix_rate, pCoord->sigma95_NEU[0], pCoord->sigma95_NEU[1], pCoord->sigma95_NEU[2]
                        , pCoord->xyz_itrf2020[0], pCoord->xyz_itrf2020[1], pCoord->xyz_itrf2020[2]
                        , pCoord->blh_itrf2020[0] * R2D, pCoord->blh_itrf2020[1] * R2D, pCoord->blh_itrf2020[2]
                        , pCoord->dN, pCoord->dE, pCoord->dU
                    );

            }
            if (fOUT) fclose(fOUT);
#endif

            plot_sol_diff(pCoords->first, pCoords->second, sN, sE, sU);
        }
    }
    return ret;
}

int main(int argc, const char *argv[])
{
    if (argc < 2)
    {
        printhelp();
    }
    else 
    {
        coord_t coord;
        std::map<std::string, coord_t> mCoords;
        if (strstr(argv[1], "update_opus_file") && coord.read_from_opus_file(argv[2]))
        {
            
        }
        else if (strstr(argv[1], "update_opus_list"))
        {
            FILE* fCSV = fopen(argv[2], "r");
            char buffer[255] = { 0 };
            while (fCSV && !feof(fCSV) && fgets(buffer, sizeof(buffer), fCSV))
            {
                char* temp = strrchr(buffer, '\n'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, '\r'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, ';'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, '#'); if (temp) temp[0] = '\0';
                if (coord.read_from_opus_file(buffer))
                {

                }
                coord.reset();
            }
            if (fCSV) fclose(fCSV);
        }
        else if (strstr(argv[1], "update_nrcan_file") && coord.read_from_nrcan_file(argv[2]))
        {
        }
        else if (strstr(argv[1], "update_nrcan_list"))
        {
            FILE* fCSV = fopen(argv[2], "r");
            char buffer[255] = { 0 };
            while (fCSV && !feof(fCSV) && fgets(buffer, sizeof(buffer), fCSV))
            {
                char* temp = strrchr(buffer, '\n'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, '\r'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, ';'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, '#'); if (temp) temp[0] = '\0';
                if (coord.read_from_nrcan_file(buffer))
                {

                }
                coord.reset();
            }
            if (fCSV) fclose(fCSV);
        }
        else if (strstr(argv[1], "update_solu_file"))
        {
            FILE* fCSV = fopen(argv[2], "r");
            char buffer[1024] = { 0 };
            while (fCSV && !feof(fCSV) && fgets(buffer, sizeof(buffer), fCSV))
            {
                char* temp = strrchr(buffer, '\n'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, '\r'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, ';'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, '#'); if (temp) temp[0] = '\0';
                if (coord.read_scsv_file(buffer))
                {
                    coord.convert_coord();
                }

                coord.reset();
            }
            if (fCSV) fclose(fCSV);
        }
        else if (strstr(argv[1], "update_esti_list"))
        {
            FILE* fCSV = fopen(argv[2], "r");
            char buffer[1024] = { 0 };
            while (fCSV && !feof(fCSV) && fgets(buffer, sizeof(buffer), fCSV))
            {
                char* temp = strrchr(buffer, '\n'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, '\r'); if (temp) temp[0] = '\0';
                if (coord.read_estimation_list(buffer))
                {
                    double rms3D = sqrt(coord.sigma95_NEU[0] * coord.sigma95_NEU[0] + coord.sigma95_NEU[1] * coord.sigma95_NEU[1] + coord.sigma95_NEU[2] * coord.sigma95_NEU[2]);
                    printf("OPUS ,%s,%14s,%s,%s,%17s,%10.5f,%14.4f,%14.4f,%14.4f,%6.2f,%10.6f,%10.6f,%10.6f,%10.6f\n"
                        , coord.code.c_str(), coord.name.c_str(), coord.stime.c_str(), coord.ctime.c_str()
                        , coord.coord_name_itrf2020.c_str(), coord.epoch_itrf2020, coord.xyz_itrf2020[0], coord.xyz_itrf2020[1], coord.xyz_itrf2020[2], coord.amb_fix_rate, coord.sigma95_NEU[0], coord.sigma95_NEU[1], coord.sigma95_NEU[2], rms3D);

                    ecef2pos_(coord.xyz_itrf2020, coord.blh_itrf2020, RE_GRS80, FE_GRS80);

                    //coord.convert_coord();

                    mCoords[coord.get_date_time()] = coord;
                }

                coord.reset();
            }
            if (fCSV) fclose(fCSV);
        }
        else if (strstr(argv[1], "update_hist_list"))
        {
            /*
            Station,Date(GNSS Time),Computing time(UTC Time),Platform,Rejected Epochs,Fixed Ambiguities,Sigmas(95%) North,Sigmas(95%) East,Sigmas(95%) Height,ITRF2020 name,epoch,x,y,z,sigma95X,sigma95Y,sigma95Z,ITRF2015 name,epoch,x,y,z,vx,vy,vz,WGS84 name,x,y,z,vx,vy,vz,Regional name,epoch,x,y,z,vx,vy,vz
            */
            FILE* fCSV = fopen(argv[2], "r");
            char buffer[1024] = { 0 };
            while (fCSV && !feof(fCSV) && fgets(buffer, sizeof(buffer), fCSV))
            {
                char* temp = strrchr(buffer, '\n'); if (temp) temp[0] = '\0';
                temp = strrchr(buffer, '\r'); if (temp) temp[0] = '\0';
                if (coord.read_history_list(buffer))
                {
                    double rms3D = sqrt(coord.sigma95_NEU[0] * coord.sigma95_NEU[0] + coord.sigma95_NEU[1] * coord.sigma95_NEU[1] + coord.sigma95_NEU[2] * coord.sigma95_NEU[2]);
                    printf("OPUS ,%s,%14s,%s,%s,%17s,%10.5f,%14.4f,%14.4f,%14.4f,%6.2f,%10.6f,%10.6f,%10.6f,%10.6f\n"
                        , coord.code.c_str(), coord.name.c_str(), coord.stime.c_str(), coord.ctime.c_str()
                        , coord.coord_name_itrf2020.c_str(), coord.epoch_itrf2020, coord.xyz_itrf2020[0], coord.xyz_itrf2020[1], coord.xyz_itrf2020[2], coord.amb_fix_rate, coord.sigma95_NEU[0], coord.sigma95_NEU[1], coord.sigma95_NEU[2], rms3D);

                    //coord.convert_coord();
                }
                coord.reset();
            }
            if (fCSV) fclose(fCSV);
        }
        else if (strstr(argv[1], "check_station_solu") && argc > 2)
        {
            std::map<std::string, std::vector<coord_t>> mCoordsAll;
            for (int i = 2; i < argc; ++i)
            {
                read_ccsv_file_all(argv[i], mCoordsAll);
            }
            process_data_shift(mCoordsAll);
        }
        else if (strstr(argv[1], "check_estimate_list") && argc > 2)
        {
            std::map<std::string, std::vector<coord_t>> mCoordsAll;
            for (int i = 2; i < argc; ++i)
            {
                read_ccsv_file1(argv[i], mCoordsAll);
            }
            process_data_shift(mCoordsAll);
        }
    }
    return 0;
}
