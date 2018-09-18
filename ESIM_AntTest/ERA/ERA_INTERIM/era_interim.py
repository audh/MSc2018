#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

server.retrieve({
 	"class": "ei",              # do not change
    "dataset": "interim",       # do not change
    "expver": "1",              # do not change
    "stream": "oper",           # monthly means of daily means
    "type": "an",               # analysis (versus forecast, fc)
    "levtype": "sfc",           
    "param": "500034.128",         # here: sea surface temperature (param 34) and mean sea level pressure (param 151)
    "time": "0000/0600/1200/1800",
	"grid": "0.75/0.75",        # horizontal resolution of output in degrees lat/lon
	"area": "-55/25/-75/35",    
	"format": "netcdf",         # get output in netcdf; only works with regular grids; for GRIB remove this line
    "date": "2017-01-01/to/2017-01-31",       # dates, set automatically from above
    "target": "qa_test"           # output file name, set automatically from above
})
 
#/500035/167.128/164.128/260087/165.128/166.128/3059
