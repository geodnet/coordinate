# Coordinate  
Coordinate is a utility used by GEODNET to extract NOAA/NGS OPUS && NRCAN PPP ITRF solution to get ITRF2020, ITRF2014, Regional (NAD83,ETRF2010, etc.) and WGS84 coordinate and velocity (based on ITRF2020 velocity field). 

# Coordinate Estimation Procedure of GEODNET Base Station
## New Station
If a new station is online for 3 continuous hours, then we will collect 30s sample interval RINEX data and call NRCAN PPP<sup>1</sup> API. If the NRCAN PPP return to a reasonable solution (as shown as Table 1), then the station is added to the RTK service list. Otherwise, the same procedure will be repeated every hour using the last 3 hours' rinex data. For a new station, another NRCAN PPP solution will be achieved once it has a full day (in GPS time). The station coordinate of the base station will be updated based on the NRCAN PPP solution quality.

Table 1. Solution to turn a new station status from ONLINE to ACTIVE  


## Relocated Station 
If a station is online again after offline for more than 30 minutes, then the SPP (Standard Point Positioning) solution will be used to compare with the current PPP solution, if there is significant changes, it will be treated as a New Station<sup>2</sup>.  

Table 2. Solution to turn an offline station from ONLINE to ACTIVE    

## Existing station
The base station coorinate of all stations will be re-estimated every 15 days, and the station coordinate will be updated based on the NRCAN PPP solution quality. Thus the station coordinate should be always within 1 month of current epoch.  

## How to get country code for each region
Coordinate use get_code.js to get the country code (3 letter) to determine the regional coordinate system.

# How to get the ITRF2020 velocity based on ITRF2020 velocity field
Coordinate use htdp (NOAA/NGS) to compute the tectonic plates velocity in ITRF2020, then to project coordinate into various epochs. This only works for linux platform.

## Regional Geodetic Coordinate System
The regional Geodetic Coordinate System adopted by GEODNET can be found at https://github.com/geodnet/geodnet_rtk_service_coordinate_system.git  

# Installation
Before run, need to install node js and gfortran, the npm install and make htdp binary.

<sup>1</sup> Thanks for NRCAN team to allow us to use NRCAN PPP (https://webapp.csrs-scrs.nrcan-rncan.gc.ca/geod/tools-outils/ppp.php).  
<sup>2</sup> This method can only detect the station relocated to another location, small relocation of the antenna will not be detected.  Currently GEODNET is developing the procedure to detect the cm-level station movement of the base station antenna.
