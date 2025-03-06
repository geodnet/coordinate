
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


int main(int argc, const char *argv[])
{
    if (argc < 2)
    {
        printhelp();
    }
    else 
    {
        coord_t coord;
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
    }
    return 0;
}
