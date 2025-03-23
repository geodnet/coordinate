
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <map>
#include <algorithm>
#include <filesystem>

#include "coord.h"
#include "coord_sol.h"

/* help text -----------------------------------------------------------------*/
static const char* help[] = {
" coordinate ",
" coordinate update_opus_file opus_sol_file      => update OPUS solution",
" coordinate update_opus_list opus_sol_list      => update OPUS solution in a batch",
" coordinate update_nrcan_file nrcan_sol_file    => update NRCAN solution",
" coordinate update_nrcan_list nrcan_sol_list    => update NRCAN solution in a batch",
" coordinate update_solu_file fname              => update solution in a file"
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

struct rrr_t
{
    std::string name;
    std::string code;
    double token;
    int rate;
    int count;
    rrr_t()
    {
        token = 0;
        count = rate = 0;
    }
    bool operator==(const rrr_t& s)
    {
        return name == s.name;
    }
    bool operator<(const rrr_t& s)
    {
        return name != s.name ? rate < s.rate : false;
    }
    bool operator<(double r)
    {
        return rate < r;
    }
};
static int check_rrr(int argc, const char* argv[])
{
    if (argc < 2) return 0;
    char buffer[1024] = { 0 };
    char* val[100];
    FILE* fLOG = fopen(argv[1], "r");
    std::map<std::string, rrr_t> mRRR;
    double sum0 = 0;
    double sum1 = 0;
    while (fLOG && !feof(fLOG) && fgets(buffer, sizeof(buffer), fLOG))
    {
        int num = parse_fields(buffer, val, ',');
        if (num > 3)
        {
            rrr_t rrr;
            rrr.name = std::string(val[0]);
            rrr.token = atof(val[1]);
            rrr.rate = (int)(atof(val[3]) * 10000);
            sum0 += rrr.token;
            int pas_verion = atoi(val[4]);
            if (pas_verion < 2) continue;
            sum1 += rrr.token;
            mRRR[rrr.name] = rrr;
        }
    }
    if (fLOG) fclose(fLOG);
    for (int i = 2; i < argc; ++i)
    {
        fLOG = fopen(argv[i], "r");
        while (fLOG && !feof(fLOG) && fgets(buffer, sizeof(buffer), fLOG))
        {
            if (strstr(buffer, "mountpoint")) continue;
            int num = parse_fields(buffer, val, ',');
            if (num > 3)
            {
                std::string name = std::string(val[0]);
                std::string code = std::string(val[1]);
                std::map<std::string, rrr_t>::iterator pRRR = mRRR.find(name);
                if (pRRR != mRRR.end())
                {
                    pRRR->second.code = code;
                    pRRR->second.count++;
                }
            }
        }
        if (fLOG) fclose(fLOG);
    }
    std::vector< rrr_t> vRRR, vRRR1, vRRR2;
    for (std::map<std::string, rrr_t>::iterator pRRR = mRRR.begin(); pRRR != mRRR.end(); ++pRRR)
    {
        if (pRRR->second.count > 0)
        {
            vRRR.push_back(pRRR->second);
            if (pRRR->second.code == "USA")
            {
                vRRR1.push_back(pRRR->second);
            }
            else
            {
                vRRR2.push_back(pRRR->second);
            }
        }
    }
    std::sort(vRRR.begin(), vRRR.end());
    fLOG = fopen("aaa.csv", "w");
    for (std::vector<rrr_t>::iterator pRRR = vRRR.begin(); pRRR != vRRR.end(); ++pRRR)
    {
        if (fLOG) fprintf(fLOG, "%s,%s,%8.3f,%8.3f,%i\n", pRRR->name.c_str(), pRRR->code.c_str(), pRRR->rate / 100.0, pRRR->token, pRRR->count);
    }
    if (fLOG) fclose(fLOG);
    std::sort(vRRR1.begin(), vRRR1.end());
    fLOG = fopen("aaa_usa.csv", "w");
    for (std::vector<rrr_t>::iterator pRRR = vRRR1.begin(); pRRR != vRRR1.end(); ++pRRR)
    {
        if (fLOG) fprintf(fLOG, "%s,%s,%8.3f,%8.3f,%i\n", pRRR->name.c_str(), pRRR->code.c_str(), pRRR->rate / 100.0, pRRR->token, pRRR->count);
    }
    if (fLOG) fclose(fLOG);
    std::sort(vRRR2.begin(), vRRR2.end());
    fLOG = fopen("aaa_oth.csv", "w");
    for (std::vector<rrr_t>::iterator pRRR = vRRR2.begin(); pRRR != vRRR2.end(); ++pRRR)
    {
        if (fLOG) fprintf(fLOG, "%s,%s,%8.3f,%8.3f,%i\n", pRRR->name.c_str(), pRRR->code.c_str(), pRRR->rate / 100.0, pRRR->token, pRRR->count);
    }
    if (fLOG) fclose(fLOG);
    std::vector<rrr_t>::iterator pRRR_999 = std::lower_bound(vRRR.begin(), vRRR.end(), 9990);
    std::vector<rrr_t>::iterator pRRR_995 = std::lower_bound(vRRR.begin(), vRRR.end(), 9950);
    std::vector<rrr_t>::iterator pRRR_990 = std::lower_bound(vRRR.begin(), vRRR.end(), 9900);
    std::vector<rrr_t>::iterator pRRR_950 = std::lower_bound(vRRR.begin(), vRRR.end(), 9500);
    std::vector<rrr_t>::iterator pRRR_100 = std::lower_bound(vRRR.begin(), vRRR.end(), 10000);
    int count = (int)vRRR.size();
    int loc999 = (int)(pRRR_999 - vRRR.begin());
    int loc995 = (int)(pRRR_995 - vRRR.begin());
    int loc990 = (int)(pRRR_990 - vRRR.begin());
    int loc950 = (int)(pRRR_950 - vRRR.begin());
    int loc100 = (int)(pRRR_100 - vRRR.begin());

    printf("ALL,%6i,%6i,%6i,%6i,%6i,%6i,%7.2f(95.0),%7.2f(99.0),%7.2f(99.5),%7.2f(99.9),%7.2f(100)\n", count, count - loc950, count - loc990, count - loc995, count - loc999, count - loc100, 100.0 - loc950 * 100.0 / count, 100.0 - loc990 * 100.0 / count, 100.0 - loc995 * 100.0 / count, 100.0 - loc999 * 100.0 / count, 100.0 - loc100 * 100.0 / count);

    pRRR_999 = std::lower_bound(vRRR1.begin(), vRRR1.end(), 9990);
    pRRR_995 = std::lower_bound(vRRR1.begin(), vRRR1.end(), 9950);
    pRRR_990 = std::lower_bound(vRRR1.begin(), vRRR1.end(), 9900);
    pRRR_950 = std::lower_bound(vRRR1.begin(), vRRR1.end(), 9500);
    pRRR_100 = std::lower_bound(vRRR1.begin(), vRRR1.end(), 10000);
    count = (int)vRRR1.size();
    loc999 = (int)(pRRR_999 - vRRR1.begin()), loc995 = (int)(pRRR_995 - vRRR1.begin()), loc990 = (int)(pRRR_990 - vRRR1.begin()), loc950 = (int)(pRRR_950 - vRRR1.begin()), loc100 = (int)(pRRR_100 - vRRR1.begin());

    printf("USA,%6i,%6i,%6i,%6i,%6i,%6i,%7.2f(95.0),%7.2f(99.0),%7.2f(99.5),%7.2f(99.9),%7.2f(100)\n", count, count - loc950, count - loc990, count - loc995, count - loc999, count - loc100, 100.0 - loc950 * 100.0 / count, 100.0 - loc990 * 100.0 / count, 100.0 - loc995 * 100.0 / count, 100.0 - loc999 * 100.0 / count, 100.0 - loc100 * 100.0 / count);

    pRRR_999 = std::lower_bound(vRRR2.begin(), vRRR2.end(), 9990);
    pRRR_995 = std::lower_bound(vRRR2.begin(), vRRR2.end(), 9950);
    pRRR_990 = std::lower_bound(vRRR2.begin(), vRRR2.end(), 9900);
    pRRR_950 = std::lower_bound(vRRR2.begin(), vRRR2.end(), 9500);
    pRRR_100 = std::lower_bound(vRRR2.begin(), vRRR2.end(), 10000);
    count = (int)vRRR2.size();
    loc999 = (int)(pRRR_999 - vRRR2.begin()), loc995 = (int)(pRRR_995 - vRRR2.begin()), loc990 = (int)(pRRR_990 - vRRR2.begin()), loc950 = (int)(pRRR_950 - vRRR2.begin()), loc100 = (int)(pRRR_100 - vRRR2.begin());

    printf("EU ,%6i,%6i,%6i,%6i,%6i,%6i,%7.2f(95.0),%7.2f(99.0),%7.2f(99.5),%7.2f(99.9),%7.2f(100)\n", count, count - loc950, count - loc990, count - loc995, count - loc999, count - loc100, 100.0 - loc950 * 100.0 / count, 100.0 - loc990 * 100.0 / count, 100.0 - loc995 * 100.0 / count, 100.0 - loc999 * 100.0 / count, 100.0 - loc100 * 100.0 / count);

    return 0;
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
    }
    return 0;
}
