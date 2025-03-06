Coordinate
Coordinate is a utility to extract NOAA/NGS OPUS && NRCAN PPP ITRF solution to get ITRF2020, ITRF2014, Regional (NAD83,ETRF2010, etc.) and WGS84 coordinate and velicity.

Coordinate use get_code.js to get the country code (3 letter) to determine the regional coordinate system.

Coordinate use htdp (NOAA/NGS) to compute the tectonic plates velocity in ITRF2020, then to project coordinate into various epochs. This only works for linux platform.

Before run, need to install node js and gfortran, the npm install and make htdp binary.


