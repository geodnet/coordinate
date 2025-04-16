# Coordinate  
Coordinate is a utility used by GEODNET to extract NOAA/NGS OPUS && NRCAN PPP ITRF solution to get ITRF2020, ITRF2014, Regional (NAD83,ETRF2010, etc.) and WGS84 coordinate and velocity (based on ITRF2020 velocity field). 

# Coordinate Estimation Procedure of GEODNET Base Station
## New Station
When a new station is online and it will be marked as ONLINE as its status. After 3 hours' continuous online, the 30s sample interval RINEX data will be uploaded to NRCAN PPP<sup>1</sup> for precise point positioning. If the NRCAN PPP return to a reasonable solution (***3D Sigma in 95% < 0.05 m***), then the station is added to the RTK service list, the station status will change from ONLINE to ACTIVE. Otherwise, the same procedure will be repeated every hour using the last 3 hours' rinex data. For a new station, another NRCAN PPP solution will be achieved once it has a full day (in GPS time). The station coordinate of the base station will be updated based on the NRCAN PPP solution quality.

## Offline Station  
If a station goes offline for ***more than 30 minutes***, it will be marked as offline

## Relocated Station 
If a station is online again after offline for more than ***30 minutes***, it will be treated as a New Station, and NRCAN PPP will be used to get the new solution. If the new NRCAN PPP solution meet the requirement (***3D Sigma in 95% < 0.05 m***), the station status will switch from ONLINE to ACTIVE. Then an coordinate update strategy will be applied to determine which solution (new NRCAN PPP solution vs existing NRCAN PPP solution) will be used.

Table 1. Coordinate Update strategy for possible relocated stations
|#|Condition1|Condition2|Solution Results|Station Relocated or not|
|:---:|:---:|:---:|:---:|:---:|
|1|3D sigma in 95% of existing solution >= 3D sigma in 95% of new solution|NONE|Use the new solution|NONE|
|2|3D sigma in 95% of existing solution <  3D sigma in 95% of new solution|3D coordinate difference between existing and new solution > ( SQRT( SQ(3D sigma in 95% of existing solution) + SQ(3D sigma in 95% of new solution) )|Use the new solution|YES|
|3|3D sigma in 95% of existing solution <  3D sigma in 95% of new solution|3D coordinate difference between existing and new solution <=( SQRT( SQ(3D sigma in 95% of existing solution) + SQ(3D sigma in 95% of new solution) )|Use existing solution|NO|

It should be noted that SQRT = **SQ**uare **R**oo**T**, SQ = **SQ**uare.   

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
