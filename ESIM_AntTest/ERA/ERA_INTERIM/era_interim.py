#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
    
server = ECMWFDataServer()
    
server.retrieve({
    'stream'    : "oper",
    'levtype'   : "sfc",
    'param'     : "228.128",
    'dataset'   : "interim",
    'step'      : "0",
    'grid'      : "0.75/0.75",
    'time'      : "00/06/12/18",
    'date'      : "2017-12-01/to/2017-12-31",
    'type'      : "an",
    'class'     : "ei",
    'area'      : "-55/25/-75/35",
    'format'    : "netcdf",
    'target'    : "ppt.nc"
 })
